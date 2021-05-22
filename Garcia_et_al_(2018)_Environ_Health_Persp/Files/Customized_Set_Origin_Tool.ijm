//----------------------- Set Origin Tool (customized) ---------------------- 
//
// Based on the 'Set Origin Tool' 
// from the following ImageJ thread:
// http://imagej.1557.x6.nabble.com/How-to-SET-a-point-as-the-new-ZERO-AXIS-POINT-td4651547.html
// 
// Customized by Cheryl L. Dunham
// Budding Bioinformatician
// Tanguay Lab | Sinnhuber Aquatic Research Laboratory (SARL)
// Oregon State University
// Corvallis, Oregon
// phone: (541) 737-3608
// email: dunhamcg@gmail.com
//
// Modified to work with the modified "Click Coordinates Tool". 
// This macro will set the origin of an image on click. The point will
// be set as (0,0) on an (x,y) pixel-based coordinate system.

  macro "Set Origin Tool - C00fL808fL08f8" {;
    // added below line to ensure scale is set in pixels  
     run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel"); 
     getCursorLoc(x, y, z, flags); 
     run("Properties...", "origin="+ x+","+y); 
     showStatus("Origin set to "+x+","+y); 
   } 
