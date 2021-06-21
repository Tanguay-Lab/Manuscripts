//------------------- Customized Click Coordinates Tool ------------------- 
//
// Based on the 'C l i c k   C o o r d i n a t e s   T o o l' 
// from the ImageJ Developer Resources Macros page:
// https://imagej.nih.gov/ij/macros/tools/ClickCoordinatesTool.txt
// 
// Customized by Cheryl L. Dunham
// Budding Bioinformatician
// Tanguay Lab | Sinnhuber Aquatic Research Laboratory (SARL)
// Environmental and Molecular Biology
// Oregon State University
// Corvallis, Oregon
// phone: (541) 737-3608
// email: dunhamcg@gmail.com
//
// On each click into an image, the x, y coordinates of the point are 
// written into the "Results" window. The point can be also marked
// in the image (destructively). This tool can handle scaled images 
// (also with nontrivial pixel aspect ratio). Double click on the tool
// icon to display the options dialog box. The "Invert Y" option in 
// Analyze > Set Measurements is supported.
// 
// Modified to work with the "Customized Set Origin Tool". 
// https://github.com/Tanguay-Lab/Manuscripts/blob/main/Garcia_et_al_(2018)_Environ_Health_Persp/Files/Customized_Set_Origin_Tool.ijm


  // removed  var outputScaled = 1;     report raw coordinates (pixels) if false
  var drawPoints = 0;      // draw cross at position of click
  var drawNumbers = 0;     // draw line number for each click

  macro 'Click Coordinates Tool - C000P515335150P5a595775950D46D64P88ab0D8bDa8Pe8cc0Pc8c90D9fDbfDdf' {
     requires("1.37e");
     getCursorLoc(x, y, z, flags); // retrieve pixel coordinates on click
     if (drawPoints || drawNumbers) setupUndo();
     if (drawPoints) {
        setLineWidth(1);
        tickLength = 3;	// the "radius" of the crosses marking the points
        drawLine(maxOf(x-tickLength,0),y, minOf(x+tickLength,getWidth()-1), y);
        drawLine(x,maxOf(y-tickLength,0), x, minOf(y+tickLength,getHeight()-1));
     }
     if (drawNumbers) {
        setFont("SansSerif",9);
        if (drawPoints) {
           setJustification("left");
           xText = x + tickLength + 1;
        } else {
           setJustification("center");
           xText = x + 1;
        }
        drawString(nResults+1, xText, y+6);
     }
     invertY = parseInt(call("ij.plugin.filter.Analyzer.getMeasurements"))&4096!=0;
     if (invertY) y = getHeight() - y - 1;
     xScale = 1;
     yScale = 1;
	 
	 // removed outputScaled if else statement
	 // replaced with below three lines of code
	 toScaled(scaledx, scaledy);  // gets "Set Origin Tool"-defined origin offset
	 newx = x + scaledx; // sets the x coordinate in relation to the user-defined origin
	 newy = -(y + scaledy); // sets the y coordinate in relation to the user-defined origin
	 
	 // below 2 lines of code replace the original setResult() lines
	 setResult("X", nResults, newx); 
     setResult("Y", nResults-1, newy);
	 updateResults();

  }

  macro 'Click Coordinates Tool Options...' {
     requires("1.37e");
     Dialog.create("Click Coordinates Tool Options");
     Dialog.addCheckbox("Draw Cross at Each Clicked Point", drawPoints);
     Dialog.addCheckbox("Write Point Number at Each Clicked Point", drawNumbers);
     Dialog.show();
	 // removed outputScaled = Dialog.getCheckbox();
     drawPoints = Dialog.getCheckbox();
     drawNumbers = Dialog.getCheckbox();
  }
