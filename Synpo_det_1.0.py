from ij import IJ, WindowManager, ImagePlus
from ij.gui import Roi
from ij.plugin.frame import RoiManager
from ij.process import ImageProcessor
from ij.measure import ResultsTable, Measurements
from ij.plugin.filter import ParticleAnalyzer
from java.lang import Double
from ij import IJ
from ij.gui import Roi, PolygonRoi
from ij.measure import ResultsTable
from ij.io import DirectoryChooser
import os
from ij.io import RoiDecoder
from javax.swing import JOptionPane
from ij.text import TextWindow
import os
from ij.gui import ShapeRoi
from ij.process import ByteProcessor
from ij.process import Blitter

### user input ###
na_im_folder="images"   ##type in the folder name containing your image
na_ROI_folder="ROI_synpodet"    ##type in the folder name containing your image
## folder selection for data and results storage
dc=DirectoryChooser("select folder contains your images and ROIs")
rsdc=DirectoryChooser("select folder to save the detected ROIs")
resdc=DirectoryChooser("select folder to save the detected results table")
pa_log=DirectoryChooser("select folder to save log file")

## threshold for partical detection (set lower threshold )
pa_thre=0.6
area_thr=0.5
                                                           
### main code ######################################################################
### do not change!!!################################################################
## define functions to ask user if all windows should be closed ##
def close_all_windows():
    # Close all image windows
    #IJ.run("Close All")
    WindowManager.closeAllWindows()
    # Close all results tables
    non_image_wins=WindowManager.getNonImageTitles()
    
    for t in non_image_wins:
        if not (t.startswith("Recorder") or t.startswith("Untitled") or "Script Editor" in t):
            win= WindowManager.getWindow(t)
            if isinstance(win, TextWindow):
                win.close()  
            
## define a function to ask if user want to proceed or not.
def askYesNo(title, message):

    # Show a YES/NO dialog
    response = JOptionPane.showConfirmDialog(None, message, title, JOptionPane.YES_NO_OPTION)
    
    # Interpret the user's response
    if response == JOptionPane.YES_OPTION:
        return True
    elif response == JOptionPane.NO_OPTION or response == JOptionPane.CLOSED_OPTION:
        return False 


##generate directory lists                                                        
main_directory = dc.getDirectory()
roi_save = rsdc.getDirectory()
results_save= resdc.getDirectory()
log_folder=pa_log.getDirectory()

## define a function to process batch of images
def batch_analysis_synpo():
    IJ.log("=== Analysis start ===")
    IJ.log("intensity threshold for detection: "+str(pa_thre))
    IJ.log("area threshold for detection: "+str(area_thr))
    IJ.log("Initialising parameters for image handling...")
    image = None
    mask_image= None
    ip_edge = None
    num_pro_images=0
    IJ.log("Sorting directories...")
    IJ.log("your main directory is: " + main_directory) ## log which directory is selected
    ## check if the directory exist
    if not main_directory:   
        IJ.showMessage("No directory was chosen")
        return
    ## creat path to find images and ROIs that need to be processed 
    image_directory = os.path.join(main_directory, na_im_folder)
    roi_directory = os.path.join(main_directory, na_ROI_folder)
    ## log directories and start processing
    IJ.log("ROI directory is: " + str(roi_directory))
    IJ.log("image directory is: " + str(image_directory))
    IJ.log("Results saving dir: "+ str(results_save))
    IJ.log("Detected ROIs saving dir: "+ str(roi_save))
    IJ.log("###start processing...###")
    
    roi_manager = RoiManager.getInstance() ## check status of ROI manager. If there is no ROI manager, open a new one                              
    if roi_manager is None:
        IJ.log("ROI manager (rm) is closed, no ROIs in rm. Opening new rm ")
        roi_manager = RoiManager()
    
    ## list all images in the image_directory, creat name for saving detected ROIs in zip file and corresponding results in csv format ##    
    for f in os.listdir(image_directory):
        if f.endswith(".tif"):
            IJ.log("TIF file detected!")
            num_pro_images+=1
            ##creating path to save detected ROIs and corresponding results
            roi_save_path = os.path.join(roi_save, f + ".zip")
            results_save_path = os.path.join(results_save, f + "_Results.csv")

            ## Note!!! the name of the ROI file should be the same as images
            roi_filename = os.path.splitext(f)[0] + ".roi"                                                                
            #zip_filename = os.path.splitext(f)[0] + ".zip"   
            ## log ROI names
            #IJ.log("ROI name (roi): " + str(roi_filename))
            #IJ.log("ROI name (zip): " + str(zip_filename))
            IJ.log("opening image...")
            ## open image if it is in tif format ##
            image = IJ.openImage(os.path.join(image_directory, f))  ## open image
            IJ.log("image being processed: " + str(f)) ## log image that is being processed
            image.show()  ## show image
            
            #IJ.log("!!Diagnose_max int: " + str(int_max)) ## self diagnosis of maximum intensity on opened image
            
            ## open corresponding ROI and load into rm ##
            if not os.path.exists(os.path.join(roi_directory, roi_filename)):  ## check if ROI directory exist ##
                IJ.showMessage("Error!! (no ROI directory found)")
                return
            
            IJ.log("matching ROI found!")  ## return matching ROI found in the log file
            IJ.log("Loading ROI file: " + str(roi_filename))  ## log name of the ROI being loaded    
            roi_manager.runCommand("Open", os.path.join(roi_directory, roi_filename)) ## open the ROI and load into rm ##
            num_rois = roi_manager.getCount() ## count how many ROIs in rm
            IJ.log("Number of ROIs loaded in rm: " + str(num_rois)) ## log numer of ROIs ##
            ## only if there are ROIs in rm, start image processing ##
            if not num_rois > 0:
                IJ.showMessage("rm is empty")
                return
            image = IJ.getImage() ## select front image 
            IJ.log("name of image for creating mask: "+ image.getTitle())
            roi = roi_manager.getRoi(0) ## get ROI from rm
            ## check if roi is polygon shape
            if roi.getType() != Roi.POLYGON:
                raise TypeError("The ROI is not a polygon")
            ## Get the polygon points of roi
            xpoints = roi.getPolygon().xpoints
            ypoints = roi.getPolygon().ypoints
            npoints = roi.getPolygon().npoints
            
            image.setRoi(roi) ## set ROI on image
            cropped_image=image.crop()  ##crop image using selected ROI
            
            ## Create a new polygon ROI for the cropped image
            xpoints_translated = [x - roi.getBounds().x for x in xpoints]
            ypoints_translated = [y - roi.getBounds().y for y in ypoints]

            new_roi = PolygonRoi(xpoints_translated, ypoints_translated, npoints, Roi.POLYGON)
            cropped_image.setRoi(new_roi)
            roi_manager.addRoi(new_roi)
            cropped_image.show()
            
            ip = cropped_image.getProcessor() ## set cropped image to processer as ip                     
            ## Apply a threshold to generate a binary mask of the cropped image
            int_max = ip.getMax() ## measure max intensity and threshold
            ip.setThreshold(pa_thre*int_max, int_max, ImageProcessor.NO_LUT_UPDATE) ## apply threshold 
            IJ.log("max intensity: "+str(int_max)) ## log max intensity 
            IJ.log("Detection threshold (min): "+str(pa_thre*int_max)) ## log threshold applied on the image
            binary_mask = ip.createMask()  ## make a binary mask based on the applied threshold
            mask_image = ImagePlus("Binary Mask", binary_mask)  ## convert binary_mask to an actual image
            IJ.run(mask_image, "Convert to Mask", "")
            mask_image.show() ## show binary mask image
            ## Now analyze this thresholded mask_image with the Particle Analyzer
            IJ.run(mask_image, "Find Edges", "")  ## find intensity peaks from clusters or punctas
            #ip_edge = WindowManager.getCurrentImage().getProcessor()  ## Getting the mask_image after finding edges for ROI detection
            ip_edge = mask_image.duplicate()
            ## if no intensity peaks were found, show error message
            if ip_edge is None:
                raise Exception("The edge detection did not work as expected. 'ip_edge' is None.")
            ## get the new roi and check ##
            roi_c=roi_manager.getRoi(1)
            if roi_c is None:
                raise ValueError("ROI not found in ROI Manager")       
            
            ## creat blank processor ## 
            blank_processor = ByteProcessor(mask_image.getWidth(), mask_image.getHeight())
            blank = ImagePlus("Blank", blank_processor)

            bounds = roi_c.getBounds()
            roi_mask = roi_c.getMask()

            ## Get the bounding rectangle of the ROI
            blank.getProcessor().insert(roi_mask, bounds.x, bounds.y)
            ## Apply this mask to the original mask image after find edges ##
            ipb = ip_edge.getProcessor()
            ipb.copyBits(blank.getProcessor(), 0, 0, Blitter.AND)
            
            ip_edge.setProcessor(ipb)
            ip_edge.updateAndDraw()
            ip_edge.show()                   
                                                    
            ## Analyze particles to get the outlines as ROIs using polygon shape
            rt = ResultsTable() ## set results table to store measured information
            roi_manager.reset() ## clear rm 
            ParticleAnalyzer.setResultsTable(rt) ## prepare rt
            ParticleAnalyzer.setRoiManager(roi_manager) ## prepare rm
            #pa = ParticleAnalyzer(ParticleAnalyzer.ADD_TO_MANAGER, 0, rt, 0, Double.POSITIVE_INFINITY, 0.0, 1.0) ## set partical analyser for cluster and puncta detection, 
                                                                                                                         ## use polygon shape to outline the clusters based on the applied threshold
                                                                                                                         ## add detected ROIs into manager
            pa = ParticleAnalyzer(ParticleAnalyzer.ADD_TO_MANAGER | ParticleAnalyzer.SHOW_MASKS, 0, rt, 0, float('inf'), 0.0, 1.0)
            pa.analyze(ip_edge)
            rt.show("particle analyze")  ## Show results table which has all measurements

            ## Filter out small ROIs based on area size 
            roi_areas = []  ## initialise an empty list to store area size
            nu_det_rois=roi_manager.getCount() ## count number of detected clusters
            IJ.log("Number of detected ROIs is:"+str(nu_det_rois))  ## log the detected cluster number
                    ## loop through the detected ROIs and measure area size
            for index in range(roi_manager.getCount()):
                roi_d = roi_manager.getRoi(index)
                image.setRoi(roi_d)
                stats = image.getStatistics(Measurements.AREA)
                roi_areas.append(stats.area)
    
            ## Calculate average area and threshold
            IJ.log("Filtering ROIs based on area size...")
            average_area = sum(roi_areas) / len(roi_areas) if roi_areas else 0   ## calculate average area size of detected cluster
            threshold_area = area_thr * average_area           ## set threshold for area size to remove false clusters
            IJ.log("average area size of all detected ROIs: "+ str(average_area))  ## log area size
            IJ.log("area threshold: "+ str(threshold_area))  ## log threshold
            ## Identify and remove small ROIs
            rois_to_remove = [index for index, area in enumerate(roi_areas) if area < threshold_area]
            for index in sorted(rois_to_remove, reverse=True):
                roi_manager.select(index)
                roi_manager.runCommand("Delete")
                    
            ## measure clusters representing actual synpo clusters
            num_detrois=roi_manager.getCount()
            IJ.log("Number of detected ROIs after filtering: "+str(num_detrois))
            IJ.log("Measureing filtered ROIs...")
            roi_manager.runCommand("Measure")
            filtered_rt = ResultsTable.getResultsTable()
            filtered_rt.show("Filtered Results")
            roi_manager.runCommand("Save", roi_save_path)
            IJ.log("ROIs saved to: " + roi_save_path)
            filtered_rt.saveAs(results_save_path)
            IJ.log("Filtered Results Table saved to: " + results_save_path)
            IJ.log("==============================")
            IJ.log("current image finished!")
            IJ.log("clearing ROIs, results tables...")
            roi_manager.runCommand("Reset")
            rt.reset()
            filtered_rt.reset()
            IJ.log("closing images...")
            image.changes = False
            image.close()
            mask_image.changes = False
            mask_image.close()
            ip_edge.changes= False
            ip_edge.close()
            cropped_image.changes = False
            cropped_image.close()
            ## reseting all variables ##
            IJ.log("Reseting all variables...")
            int_max=0
            IJ.log("int_max after reset: "+str(int_max))
            if not int_max == 0:
                IJ.showMessage("Error! int_max is incorrect")
                return
            ip = None
            image = None
            mask_image= None
            ip_edge = None
            cropped_image = None
            average_area = 0
            threshold_area = 0
            rois_to_remove= []
            num_rois = 0
            
            IJ.log("Done!")
            IJ.log("Moving to next image...")
            IJ.log("###############################")
                         
        else:
            IJ.showMessage("Error. No tif image detected.")
            IJ.log("Error. No tif image detected.")
 
    IJ.log("===========================")
    IJ.log("Analysis finished")
    IJ.log("Number of images processed: "+str(num_pro_images))
 
    
## check if the function defination is working. 
## if so, start processing
    
if __name__ == '__main__':
    batch_analysis_synpo()
    ## ask if user want to save log file and close all windows in FIJI ##
       
    userResponse = askYesNo("Saving log file_closing all windows", "Do you want to save the log file and CLOSE all windows?")
    if userResponse:
        path_log=os.path.join(log_folder, "ImageJ_Log.txt")
        IJ.log("Saving log file to: "+str(path_log))
        log = IJ.getLog()
        with open(path_log, "w") as logfile:
            logfile.write(log)
        close_all_windows()
    else:
        IJ.showMessage("Analysis finished (log file unsaved!)")
        IJ.log("log file unsaved.")     

### this is the end! ###