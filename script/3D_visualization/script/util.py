import tifffile as tif
import numpy as np
import pandas as pd
import re
import os
import string
from pathlib import Path, PurePath
from PIL import Image, ImageSequence, UnidentifiedImageError


def pil_imopen(fname, metadata=False):
    im = Image.open(fname)

    if metadata:
        return im, pil_getmetadata(im)
    else:
        return im


def pil_imread(
    fname,
    metadata=False,
    swapaxes=False,
    ensure_4d=True,
    backup=tif.imread,
    **kwargs
):
    md = None

    import warnings
    warnings.simplefilter('ignore', UserWarning)

    try:
        im = pil_imopen(fname)
        md = pil_getmetadata(im)
        imarr = pil_frames_to_ndarray(im)
    except (ValueError, UnidentifiedImageError) as e:
        if callable(backup):
            imarr = backup(fname, **kwargs)
        else:
            raise e

    if ensure_4d and imarr.ndim == 3:
        # assumes 1 Z
        imarr = imarr[:, None, :]

    if swapaxes and imarr.ndim == 4:
        imarr = imarr.swapaxes(0, 1)

    if metadata and md:
        return imarr, md
    else:
        return imarr


def pil_getmetadata(im, relevant_keys=None):
    """
    pil_getmetadata
    ---------------
    Given a PIL image sequence im, retrieve the metadata associated
    with each frame in the sequence. Only keep metadata keys specified
    in `relevant_keys` - which will default to ones that we need such as
    channel, slice information. There are many metadata keys which are
    useless / do not change frame to frame.
    Returns: List of dicts in order of frame index.
    """

    if str(relevant_keys).lower() == 'all':
        relevant_keys = None

    elif not isinstance(relevant_keys, list):

        relevant_keys = [
            'Andor sCMOS Camera-Exposure',  # Exposure time (ms)
            'Channel',                      # Channel name (wavelength)
            'ChannelIndex',                 # Channel index (number)
            'Frame',                        # Time slice (usually not used)
            'FrameIndex',                   # Time slice index (usually not used)
            'PixelSizeUm',                  # XY pixel size in microns
            'Position',                     # Position 
            'PositionIndex',                # Position index (MMStack_PosX)
            'PositionName',                 # Position name
            'Slice',                        # Z slice
            'SliceIndex'                    # Z slice index (same as Slice)
        ]

    frame_metadata = []

    for frame in ImageSequence.Iterator(im):

        # The JSON string is stored in a key named "unknown",
        # probably because it doesn't correspond to a standard
        # TIF tag number.
        if 'unknown' in frame.tag.named().keys():
            jsstr = frame.tag.named()['unknown'][0]
            jsdict = json.loads(jsstr)

            if relevant_keys:
                # Only keep the relevant keys
                rel_dict = {
                    k: jsdict.get(k)
                    for k in relevant_keys
                }
            else:
                rel_dict = jsdict

            frame_metadata.append(rel_dict)

    return frame_metadata


def pil2numpy(im, dtype=np.uint16):

    return np.frombuffer(im.tobytes(), dtype=dtype).reshape(im.size)


def pil_frames_to_ndarray(im, dtype=np.uint16):
    """
    pil_frames_to_ndarray
    -----------------
    Given a PIL image sequence, return a Numpy array that is correctly
    ordered and shaped as (n_channels, n_slices, ...) so that we can 
    process it in a consistent way.
    To do this, we look at the ChannelIndex and SliceIndex of each frame
    in the stack, and insert them one by one into the correct position
    of a 4D numpy array.
    """
    metadata = pil_getmetadata(im)

    if not metadata:
        raise ValueError('Supplied image lacks metadata used for '
            'forming the correct image shape. Was the image not '
            'taken from ImageJ/MicroManager?')

    # Gives a list of ChannelIndex for each frame
    cinds = jmespath.search('[].ChannelIndex', metadata)
    # Gives a list of SliceIndex for each frame
    zinds = jmespath.search('[].SliceIndex', metadata)

    if (len(cinds) != len(zinds)
        or any([c is None for c in cinds])
        or any([z is None for z in zinds])
    ):
        raise ValueError('SuppliedImage lacks `ChannelIndex` or '
                         '`SliceIndex` metadata required to form '
                         'properly shaped numpy array. Was the image not '
                         'taken directly from ImageJ/MicroManager?')

    ncs = max(cinds) + 1
    nzs = max(zinds) + 1

    total_frames = ncs * nzs
    assert total_frames == im.n_frames, 'wrong shape'

    # Concatenate the channel and slice count to the XY shape in im.size
    new_shape = (ncs, nzs) + im.size

    # Make an empty ndarray of the proper shape and dtype
    npoutput = np.empty(new_shape, dtype=dtype)

    # Loop in a nested fashion over channel first then Z slice
    for c in range(ncs):
        for z in range(nzs):

            # Find the frame whose ChannelIndex and SliceIndex
            # match the current c and z values
            entry = jmespath.search(
                f'[?ChannelIndex==`{c}` && SliceIndex==`{z}`]', metadata)[0]

            # Find the *index* of the matching frame so that we can insert it
            ind = metadata.index(entry)

            # Select the matching frame
            im.seek(ind)

            # Copy the frame into the correct c and z position in the numpy array
            npoutput[c, z] = pil2numpy(im)

    return npoutput

