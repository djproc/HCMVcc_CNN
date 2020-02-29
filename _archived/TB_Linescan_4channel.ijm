input = getDirectory("Choose an input Directory");
output = getDirectory("Choose an output Directory");
setBatchMode(true);

list = getFileList(input);
for (k = 0; k < list.length; k++)
        {
			print(list[k]);
			allText = File.openAsString(input + "/" + list[k]);
			//print(allText)
			//tmp = split(fileName,".");
			// get file format {txt, csv}
			//posix = tmp[lengthOf(tmp)-1];
			// parse text by lines
			text = split(allText, "\n");
			
			run("Clear Results"); 

			// define array for points
			var x1points = newArray;
			var y1points = newArray;
			var x2points = newArray;
			var y2points = newArray;
			var C1url = newArray;
			var C2url = newArray;
			var C3url = newArray;
			var C4url = newArray;
			var ImageID = newArray;

			//define colums to collect data from
			hdr = split(text[0]);
			iImageID=0;
			// 1
			// 2
			iC1 = 3;
			iC2 = 4;
			iC3 = 5;
			iC4 = 6;
			// 7
			// 8
			//9-14
			iX1=15;
			iY1=16;
			iX2=17;
			iY2=18;

			// loading and parsing each line
			for (i = 1; i < (text.length); i++)
			{
			   line = split(text[i],",");
			   setOption("ExpandableArrays", true);   
			   x1points[i-1] = parseInt(line[iX1]);
			   y1points[i-1] = parseInt(line[iY1]);
			   x2points[i-1] = parseInt(line[iX2]);
			   y2points[i-1] = parseInt(line[iY2]);

			   C1url[i-1] = line[iC1];
			   C2url[i-1] = line[iC2];
			   C3url[i-1] = line[iC3];
			   C4url[i-1] = line[iC4];
			   
			   ImageID[i-1] = line[iImageID];

			   open(C1url[i-1]);
			   makeLine(x1points[i-1],y1points[i-1],x2points[i-1],y2points[i-1]);
			   Roi.setStrokeWidth(150);
			   profile = getProfile();
			   for (j=0; j<profile.length; j++)
			   		setResult(ImageID[i-1] + "_C1", j, profile[j]); 
			   updateResults;
			   close("*");
			   
			   open(C2url[i-1]);
			   makeLine(x1points[i-1],y1points[i-1],x2points[i-1],y2points[i-1]);
			   Roi.setStrokeWidth(150);
			   profile = getProfile();
			   for (j=0; j<profile.length; j++)
			   		setResult(ImageID[i-1] + "_C2", j, profile[j]); 
			   updateResults;
			   close("*");

			   open(C3url[i-1]);
			   makeLine(x1points[i-1],y1points[i-1],x2points[i-1],y2points[i-1]);
			   Roi.setStrokeWidth(150);
			   profile = getProfile();
			   for (j=0; j<profile.length; j++)
			   		setResult(ImageID[i-1] + "_C3", j, profile[j]); 
			   updateResults;
			   close("*");
			   
			   open(C4url[i-1]);
			   makeLine(x1points[i-1],y1points[i-1],x2points[i-1],y2points[i-1]);
			   Roi.setStrokeWidth(150);
			   profile = getProfile();
			   for (j=0; j<profile.length; j++)
			   		setResult(ImageID[i-1] + "_C4", j, profile[j]); 
			   updateResults;
			   close("*");
			}

			// Save as spreadsheet compatible text file //could just use home and take the data from there 
			path = output + "/profile_" + list[k];
			saveAs("Results", path);
        	
        }