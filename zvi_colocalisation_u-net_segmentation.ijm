# @String input_dir
# @int min_cell_area
# @int max_cell_area
# @String channel_combs
# @String clip_perc
# @String channel_colors
# @String coloc_colors
# @String segment_tile_size
# @String GPU_MB

// Macro performing colocalisation on czi files
// 1. Splits channels, moves raw data to a subdirectory
// 2. Detects positive cells in channels by segmentation with U-Net
// 3. Determines single- and double-positive cells using segmentation masks
// 4. Drops results to a csv file and draws outlines of segmentation on overlay image
// Provide network models and weights for each channel in subdirectory of input directory "models/mod[channel number]
// Bartek Tarkowski, LRB @ IIMCB, Aug 2020

// Function returning sorted unique array elements, by Rainer M. Engel, 2012
function ArrayUnique(array) {
	array 	= Array.sort(array);
	array 	= Array.concat(array, 999999);
	uniqueA = newArray();
	i = 0;	
   	while (i<(array.length)-1) {
		if (array[i] == array[(i)+1]) {
			//print("found: "+array[i]);			
		} else {
			uniqueA = Array.concat(uniqueA, array[i]);
		}
   		i++;
   	}
	return uniqueA;
}

// Function stitching czi series if needed and splitting into channels and saving them
function splitchChannels(input_dir, file_name){
	// create file name with tif ext
	tif_name = replace(file_name , ".czi" , ".tif");
	//	Open zvis to split channels
	run("Bio-Formats", "open=" + input_dir + "raw_data/" + file_name + 
	" open_all_series view=Hyperstack stack_order=XYCZT");
	// check whether multiple series present
	if (nImages == 1){
		// split stack to channels
		run("Stack to Images");
		// Loop thru open channels
		for (i=0 ; i < nImages() ; i++) {
			selectWindow( replace(file_name , ".czi" , "") + "-000" + (i+1) );
			// Make channel directory
			if ( ! File.isDirectory(input_dir + "/channels/ch" + (i+1) ) ) {
				File.makeDirectory(input_dir + "/channels/ch" + (i+1) ) ;
				}
			// Save channel image
			saveAs("Tiff", input_dir + "/channels/ch" + (i+1) + "/" + tif_name);
			};
		}else{
			getPixelSize(unit, pW, pH);		// get image scaling
			run("Grid/Collection stitching", "type=[Positions from file] order=[Defined by image metadata]" +
			" multi_series_file=" + input_dir + "raw_data/" + file_name + 
			" fusion_method=[Linear Blending]" +
			" regression_threshold=0.30 max/avg_displacement_threshold=2.50 absolute_displacement_threshold=3.50" +
			" increase_overlap=0 computation_parameters=[Save memory (but be slower)]" +
			" image_output=[Write to disk] output_directory=" + input_dir + "channels/");
			stitch_list = getFileList(input_dir + "channels/");
			for (sch = 0 ; sch < stitch_list.length ; sch++) {
				if(stitch_list[sch].startsWith("img_t1_z1_c") ) {
					curr_channel = substring(stitch_list[sch] , stitch_list[sch].length - 1);
					if ( ! File.isDirectory(input_dir + "/channels/ch" + curr_channel ) ) {
						File.makeDirectory(input_dir + "/channels/ch" + curr_channel ) ;
						}
					open(input_dir + "channels/" + stitch_list[sch] );
					run( "Set Scale...", "distance=" + (1.0 / pW) + " known=1 unit=" + unit.substring(0 , unit.length- 1) );
					run("Save");
					close();
					File.rename(input_dir + "channels/" + stitch_list[sch], 
								input_dir + "channels/ch" + curr_channel + "/" + tif_name);}
				}
			}
	run("Close All");
	};

// Function performing segmentation on channel files, filtering detected particle by size and saving rescaled masks
function segmentation(input_dir, curr_file){	
	// determine number of channels
	channels = ArrayUnique( split(channel_combs , ",+" ) );
	for (c = 0 ; c < channels.length ; c++){
		//make dir for masks for current channel
		if ( ! File.isDirectory(input_dir + "/masks/ch" + channels[c]) ) {
			File.makeDirectory(input_dir + "/masks/ch" + channels[c]);
			}
		// open tif tile for current channel
		open(input_dir + "channels/ch" + channels[c] + "/" + curr_file );
		// get image dims for rescaling after segmentation
		img_w = getWidth();
		img_h = getHeight();
		// define model and weight files
		model_files = getFileList(input_dir + "models/mod" + channels[c] + "/");
		for (m = 0; m < model_files.length; m++) {
			if (model_files[m].endsWith("modeldef.h5") ) { 
				model = input_dir + "models/mod" + channels[c] + "/" + model_files[m]; 
				}
			else if(model_files[m].endsWith("caffemodel.h5") ) { 
				weights = input_dir + "models/mod" + channels[c] + "/" + model_files[m]; 
				}
			};
		// run segmentation
		call('de.unifreiburg.unet.SegmentationJob.processHyperStack',
		'modelFilename=' + model +
		',Memory (MB):=' + GPU_MB +
		',weightsFilename=' + weights +
		',gpuId=GPU 0' +
		',useRemoteHost=false' +
		',processFolder=/home/bart/Desktop/' +
		',average=none,keepOriginal=false,outputScores=false,outputSoftmaxScores=false');

		// rescale segmentation to original image dims, threshold and filter particles by size
		if (isOpen(curr_file + " - 32-Bit - rescaled (xy) - normalized - score (segmentation)")) {
			selectWindow(curr_file + " - 32-Bit - rescaled (xy) - normalized - score (segmentation)");
			}else{
			selectWindow(curr_file + " - 32-Bit - normalized - score (segmentation)");
			}
		run("Size...", "width=" + img_w + " height=" + img_h + " depth=1 constrain average interpolation=Bilinear");
		run("Convert to Mask");
		run("Analyze Particles...", "size=" + min_cell_area + "-" + max_cell_area + " show=Overlay add");
		// rescale original image back
		if (isOpen(curr_file + " - 32-Bit - rescaled (xy) - normalized")) 
			{
			selectWindow(curr_file + " - 32-Bit - rescaled (xy) - normalized");
			}else{
			selectWindow(curr_file + " - 32-Bit - normalized");
			}
		run("Size...", "width=" + img_w + " height=" + img_h + " depth=1 constrain average interpolation=Bilinear");
		// ask user to correct segmentation
		roiManager("Show All");
		waitForUser("Manual segmentation correction", "Press Ok when corrected");
		roiManager("Show All");
		roiManager("Combine");
		run("Create Mask");
		roiManager("reset");
		// save corrected segmentation and close all images
		selectWindow("Mask");
		saveAs("Tiff", input_dir + "masks/ch" + (c+1) + "/" + curr_file);
		run("Close All");
		}	
	};

// Function performing colocalisations of selected channels and saving colocalisation masks to separate directory
function colocalisation(input_dir, curr_file, channel_combs){	
	// determine combinations of channels to colocalise and loop thru them
	ch_combs = split(channel_combs , ",");
	for (c = 0 ; ( c < ch_combs.length ) ; c++){
		// skip iteration if not colocalising combination (i.e. single channel)
		if (!ch_combs[c].contains("+")){continue;}
		// make dir for masks for current colocalisation
		if ( ! File.isDirectory(input_dir + "/colocs/ch" + ch_combs[c]) ) {
			File.makeDirectory(input_dir + "/colocs/ch" + ch_combs[c]);
			}
		// extract channels to colocalise and loop thru them
		curr_comb = split(ch_combs[c] , "+");
		for (cc = 0 ; cc < curr_comb.length ; cc++){
			open(input_dir + "masks/ch" + curr_comb[cc] + "/" + curr_file);
			run("Select All");
			run("Copy");
			// paste current mask to the first one in logical AND mode
			setPasteMode("AND");
			selectWindow(curr_file);
			run("Paste");
			}
		// filter particles by size
		selectWindow(curr_file);
		run("Analyze Particles...", "size=" + min_cell_area + "-" + max_cell_area + " show=Masks");
		// save filtered segmentation and close all images
		selectWindow("Mask of " + curr_file);
		run("Invert");
		saveAs("Tiff", input_dir + "colocs/ch" + ch_combs[c] + "/" + curr_file);	
		run("Close All");
	}	
};


// Make directories for raw and processed data
if ( ! File.isDirectory(input_dir + "/raw_data") ) {File.makeDirectory(input_dir + "/raw_data");}
if ( ! File.isDirectory(input_dir + "/channels") ) {File.makeDirectory(input_dir + "/channels");}
if ( ! File.isDirectory(input_dir + "/overlays") ) {File.makeDirectory(input_dir + "/overlays");}
if ( ! File.isDirectory(input_dir + "/masks") ) {File.makeDirectory(input_dir + "/masks");}
if ( ! File.isDirectory(input_dir + "/colocs") ) {File.makeDirectory(input_dir + "/colocs");}
	
// Get list of files from input directory
file_list = getFileList(input_dir);

// Main loop for creating channel and colocalisation masks using functions above
for (i = 0; i < file_list.length; i++){
	if (endsWith(file_list[i], ".czi")) {
		// Move original file to subdirectory
		File.rename(input_dir + file_list[i], input_dir + "raw_data/" + file_list[i]);
		// split and save channels
		print("stitching and splitting channels");
		splitchChannels(input_dir , file_list[i]);
		// extract tif file name for remaining functions
		curr_file = replace(file_list[i] , ".czi" , ".tif");
		// perform segmentation
		print("segmenting");
		segmentation(input_dir , curr_file);
		// perform colocalisation
		print("colocalising");
		colocalisation(input_dir, curr_file, channel_combs);
		}
	};
	

// Get list of files from first subdirectory with channel 1 of input directory
tif_list = getFileList(input_dir + "channels/ch1/");
// Get channel list
channel_list = getFileList(input_dir + "channels/");
// Get colocalisation list
coloc_list = getFileList(input_dir + "colocs/");
// Get channel color list and colocalisation color list from stiring joined by colons provided by the user
ch_color_list = split(channel_colors , ",");
coloc_color_list = split(coloc_colors , ",");
// Get channel combinations for outlines on overlay image
ch_combs = split(channel_combs , ",");


// Add colors to channels
for (ch = 0; ch < channel_list.length; ch++){	//loop thru channels
	if ( ! File.isDirectory(input_dir + "/overlays/ch" + (ch + 1) ) ) 
		{File.makeDirectory(input_dir + "/overlays/ch" + (ch + 1) );}
	for (t = 0; t < tif_list.length; t++){		//loop thru files in channel directory and open each
		open(input_dir + "channels/" + channel_list[ch] + "/" + tif_list[t]);
		}
	// make stack and enhance contrast according to stack histogam with user-defined clipping
	run("Images to Stack", "name=" + channel_list[ch] + " title=[] use");
	run("Enhance Contrast...", "saturated=" + clip_perc + " use");

	//add LUT, convert to RGB, split stack
	run(ch_color_list[ch]);
	run("RGB Color");
	run("Stack to Images");
	// loop thru all opened images and save them in channel subdirectory of overlays directory
	for (t = 0; t < nImages; t++){
		selectImage(t+1);
		saveAs("Tiff", input_dir + "overlays/" + channel_list[ch] + "/" + getTitle() );
		}
	// read dimensions of created overlay
	oW = getWidth();
	oH = getHeight();
	run("Close All");
	};

// Make csv file for results
File.saveString("name,ch_comb,count\n",  input_dir + "coloc_results.csv" )

// Drop results to csv file and make overlays with outlines of detected cells
for (t = 0; t < tif_list.length; t++){
	// loop thru channels and open them
	for (ch = 0 ; ch < channel_list.length ; ch++){
		open(input_dir + "overlays/" + channel_list[ch] + "/" + tif_list[t]);
		}
	// make stack of color channels, make projection to overlay them, save overlay, close others
	run("Images to Stack", "name=" + tif_list[t] + " title=[]");
	run("Z Project...", "projection=[Max Intensity]");
	if ( ! File.isDirectory(input_dir + "/overlays/merge/" ) ) 
		{File.makeDirectory(input_dir + "/overlays/merge/" );}
	saveAs("Tiff", input_dir + "/overlays/merge/" + tif_list[t]);
	close("\\Others");
	// loop thru colocalisations and open
	for (chc = 0 ; chc < ch_combs.length ; chc++){
		// open either colocalisation mask or single channel mask depending user-defined combination
		if ( ch_combs[chc].contains("+") ){open(input_dir + "colocs/ch" + ch_combs[chc] + "/" + tif_list[t]);}
		else {open(input_dir + "masks/ch" + ch_combs[chc] + "/" + tif_list[t]);}
		// resize canvas to fit overlay image
		run("Canvas Size...", "width=" + oW +" height=" + oH + " position=Center");
		// add outlines to ROI manager by analysing particles;
		run("Analyze Particles...", "size=" + min_cell_area + "-" + max_cell_area + " show=Overlay add");
		// Drop results to text file
		File.append(String.join( newArray(tif_list[t] , ch_combs[chc] , roiManager("count")) , "," ) , input_dir + "coloc_results.csv");
		// select overlay image, add overlays from ROI manager and clear ROI manager
		selectWindow(tif_list[t]);
		if (roiManager("count") > 0){
			run("From ROI Manager");
			roiManager("reset");
			// set overlay parameters and flatten image
			Overlay.setStrokeColor(coloc_color_list[chc]);
			Overlay.setStrokeWidth(4 - chc);
			Overlay.drawLabels(0);
			Overlay.flatten;
			// save overlay and close other files
			saveAs("Tiff", input_dir + "/overlays/" + tif_list[t]);
			}
		close("\\Others");
		}	
	run("Close All");
	}	