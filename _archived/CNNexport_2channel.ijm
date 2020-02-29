// to not display the pictures as they are goining, use this
setBatchMode(true);

// ask for a file to be imported
	
	fileName = File.openDialog("Select the file to import");
	allText = File.openAsString(fileName);
	tmp = split(fileName,".");
	// get file format {txt, csv}
	posix = tmp[lengthOf(tmp)-1];
	// parse text by lines
	text = split(allText, "\n");

//define your path, this is where files will be saved
	
	//path =  "B:\114\114_20180920_time course redux slide scan\_CP\v3_wholeslide_TB72only\_imgexport"
	dir = getDirectory("Select a location to save the files"); 
	
//clear results in case there is something you don't want in your output Results file

	run("Clear Results")


// define array for points and image URLS
	
	var x1points = newArray;
	var y1points = newArray;
	var x2points = newArray;
	var y2points = newArray;
	
	var nuc_Ori = newArray;
	var Xmidpoints= newArray;	
	var Ymidpoints = newArray;			
	
	var C1url = newArray;
	var C2url = newArray;
	//var C3url = newArray;
	//var C4url = newArray;
	
	var ImageID = newArray; 

//label you columns by number here (starting at 0, not 1)
	
	hdr = split(text[0]);

	iImageObjectID = 0;
	
	iXmid = 1;
	iYmid = 2;

	iC1url = 3;
	iC2url = 4;
	//iC3url = 5;
	//iC4url = 5;
				
// loading and parsing each line
	for (i = 2; i < (text.length); i++)
	{
	   //we split each line in the csv at commas 
	   line = split(text[i],",");
	   setOption("ExpandableArrays", true);   

	   //this loads the data from the csv into the arrays me setup before so that we can iterate over each line (i)
	   //these are specifically for the HCMV ac and nucleus, but you can change them if you need.

	   Xmidpoints[i-1] = parseInt(line[iXmid]);
	   Ymidpoints[i-1] = parseInt(line[iYmid]);
	   
	   C1url[i-1] = line[iC1url]; //note we do not use parseInt here because we want the URL as a string
	   C2url[i-1] = line[iC2url];
	   //C3url[i-1] = line[iC3url];
	   //C4url[i-1] = line[iC4url];
	   
	   ImageID[i-1] = line[iImageObjectID]; 

	   //now that the data is loaded from the csv, we can perform our operations
	   //open each image using the URL, draw a line with 150 width, perform a linescan and save the data to a column in the results tab 


	   open(C1url[i-1]);
	   rename("C1");
	   
	   open(C2url[i-1]);
	   rename("C2");
	   
	   //open(C3url[i-1]);
	   //rename("C3");
	   
	   run("Merge Channels...", "c1=C1 c2=C2 create");
	   run("Canvas Size...", "width=3248 height=3248 position=Center zero"); //expand canvas to 2048 + 600*2 for ROI on each side
	   makeRectangle(Xmidpoints[i-1]+300,Ymidpoints[i-1]+300,600,600); //note here we compensated for canvas size expansion
	   run("Duplicate...", "duplicate");	   			   
	   saveAs("Jpeg", dir + ImageID[i-1]+"_RGB.jpg");
	   close("*");
	   
	   //open(C4url[i-1]);
	   //run("Canvas Size...", "width=3548 height=3548 position=Center zero"); //expand canvas to 2048 + 750*2 for ROI on each side
	   //makeRectangle(Xmidpoints[i-1]+375,Ymidpoints[i-1]+375,750,750); //note here we compensated for canvas size expansion
	   //run("Duplicate...", "duplicate");	   			   
	   //rename("C4");
	   //save(dir+ImageID[i-1]+"_C4");
	   //close("*");
	   	   
	   
	   //merge all of these channles, cut out a 750x750 region and rotate to align every image on the axis of the linescan
	   //run("Merge Channels...", "c1=C1 c2=C2 c3=C3 c4=C4 create");
	   //makeRectangle(Xmidpoints[i-1]-375,Ymidpoints[i-1]-375,750,750);
	   //Roi.setStrokeWidth(2);
	   //setForegroundColor(255, 255, 255);
	   //run("Draw");
	   
	   
	   //run("Canvas Size...", "width=350 height=350 position=Center zero"); //expand canvas to 2048 + 750*2 for ROI on each side
	   //makeRectangle(Xmidpoints[i-1],Ymidpoints[i-1],50,50); //note here we compensated for canvas size expansion
	   //run("Duplicate...", "duplicate");
	   //save(dir+ImageID[i-1]);
	   //close("*");	
	
	}