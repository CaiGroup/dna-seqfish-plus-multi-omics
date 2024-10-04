# Comparing cell shape before and after FISH protocol 

Image stacks of live nuclei are in the `cropped_live_nuc_pore_padded folder`, and image stacks of the same nuclei after a mock seqFISH protocol are in the  `cropped_fix_its1_padded` folder. All stacks are cropped and padded so that the image stacks of each nuclues while live and after the mock seqFISH protocol have the same lateral dimensions.

Fit points in the nuclei for each image were obtained as described in our manuscript and are in the `filtered_live_pnts` and `fit_fix_its1` folders. Cell heights derived from these points are listed in `cell_heights.csv`.

`analyze_results_fig1.py.ipynb` contains shows the analysis of the data and plots figure panels sumarizing the data and showing how the heights were obtained using a representative cell.

`analyze_results_fig1_any_cell.py.ipynb` contains the same analysis as the preceding notebook, but generates figure panels illustrating the height measurement for all cells so readers, allowing readers to judge the accuracy of our height measurements for themselves.

`panel_rand_cell_final_run.py.ipynb` plots size projections of ten randomly chosen cells with their determined bottoms and tops annotated while live and after undergoing a mock seqFISH protocol annotated. The 10 randomly chosen cells are hard coded to be the same as shown in our revision letter, but readers may uncomment a code block to draw 10 different cells to plot.
