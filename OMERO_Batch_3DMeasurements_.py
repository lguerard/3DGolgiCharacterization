'''
Author: Laurent Guerard
Group: IMCF
Email: laurent.guerard@unibas.ch
Email bis: l.guerard42@gmail.com
Creation Date: Thursday, 7th March 2019 15:50:05
-----
Last Modified: Monday, 14th September 2020 9:18:48
Modified By: Laurent Guerard
-----
HISTORY:
Date (Y-M-D)    By        Comments
----------      --        ---------------------------------------------------------
2020-09-14		LG 		  Removed the convex hull measurements
2020-09-9       LG        Fixed the volume calibration and the threshold
2020-08-17      LG        Fixed the 3D_solidity name
2020-08-17      LG        Version to be published
'''

# ─── SCRIPT PARAMETERS ──────────────────────────────────────────────────────────

#@ String(label="username", description="please enter your username") USERNAME
#@ String(label="password", description="please enter your password", style="password") PASSWORD
#@ Integer(label="Dataset ID", description="OMERO dataset you wish to process") datasetid
#@ File(label="path to store results", style="directory", description="results of your script will be stored here") destination
#@ OpService ops
#@ UIService ui

# ─── IMPORTS ────────────────────────────────────────────────────────────────────

import os
import sys
import csv
import inspect
from itertools import izip

from ij import IJ, ImagePlus, ImageStack, CompositeImage, Prefs, WindowManager as wm
from ij.plugin import Duplicator

from pprint import pprint


# Bioformats imports
from loci.plugins import BF
from loci.plugins.in import ImporterOptions

# 3DSuite imports
from mcib3d.geom import Objects3DPopulation
from mcib3d.image3d import ImageInt, ImageHandler, Segment3DImage
from mcib3d.image3d.IterativeThresholding import TrackThreshold

# Omero Dependencies
from omero.gateway import Gateway
from omero.gateway import LoginCredentials
from omero.gateway import SecurityContext
from omero.gateway.exception import DSAccessException
from omero.gateway.exception import DSOutOfServiceException
from omero.gateway.facility import BrowseFacility
from omero.gateway.facility import DataManagerFacility
from omero.gateway.model import DatasetData
from omero.gateway.model import ExperimenterData
from omero.gateway.model import ProjectData
from omero.log import Logger
from omero.log import SimpleLogger
from omero.model import Pixels

from ome.formats.importer import ImportConfig
from ome.formats.importer import OMEROWrapper
from ome.formats.importer import ImportLibrary
from ome.formats.importer import ImportCandidates
from ome.formats.importer.cli import ErrorHandler
from ome.formats.importer.cli import LoggingImportMonitor

# java imports
from java.lang import Long
from java.lang import String
from java.lang.Long import longValue
from java.util import ArrayList
from jarray import array
from java.lang.reflect import Array
import java

# ─── FUNCTIONS ──────────────────────────────────────────────────────────────────

def openImagePlus(HOST, USERNAME, PASSWORD, groupId, imageId):
    """Open an ImagePlus from an OMERO server

    Parameters
    ----------
    HOST : str
        Adress of your OMERO server
    USERNAME : str
        Username to use in OMERO
    PASSWORD : str
        Password
    groupId : double
        OMERO group ID 
    imageId : int
        ID of the image to open
    """
    stackview = "viewhyperstack=false stackorder=XYCZT "
    datasetorg = "groupfiles=false swapdimensions=false openallseries=false concatenate=false stitchtiles=false"
    coloropt = "colormode=Default autoscale=true"
    metadataview = "showmetadata=false showomexml=false showrois=true setroismode=roimanager"
    memorymanage = "virtual=false specifyranges=false setcrop=false"
    split = " splitchannels=false splitfocalplanes=false splittimepoints=false"
    other = "windowless=true"
    options = ("location=[OMERO] open=[omero:server=%s\nuser=%s\npass=%s\ngroupID=%s\niid=%s] %s %s %s %s %s %s %s " %
               (HOST, USERNAME, PASSWORD, groupId, imageId, stackview, datasetorg, coloropt, metadataview, memorymanage, split, other))
    IJ.runPlugIn("loci.plugins.LociImporter", options)


def omeroConnect():
    """Connect to OMERO using the credentials entered

    Returns
    -------
    gateway
        OMERO gateway
    """
    # Omero Connect with credentials and simpleLogger
    cred = LoginCredentials()
    cred.getServer().setHostname(HOST)
    cred.getServer().setPort(PORT)
    cred.getUser().setUsername(USERNAME.strip())
    cred.getUser().setPassword(PASSWORD.strip())
    simpleLogger = SimpleLogger()
    gateway = Gateway(simpleLogger)
    gateway.connect(cred)
    return gateway

# List all ImageId's under a Project/Dataset


def getImageIds(gateway, datasetId):
    """List all ImageIds under a Project/Dataset

    Parameters
    ----------
    gateway : gateway
        Gateway to the OMERO server
    datasetId : int
        ID of the dataset in OMERO

    Returns
    -------
    int[]
        List of all the ImageIDs from the Dataset
    """
    browse = gateway.getFacility(BrowseFacility)
    user = gateway.getLoggedInUser()
    ctx = SecurityContext(user.getGroupId())
    ids = ArrayList(1)
    val = Long(datasetId)
    ids.add(val)
    images = browse.getImagesForDatasets(ctx, ids)
    j = images.iterator()
    imageIds = []
    while j.hasNext():
        image = j.next()
        imageIds.append(String.valueOf(image.getId()))
    return imageIds


def BFImport(indivFile):
    """Import using BioFormats

    Parameters
    ----------
    indivFile : {str}
        Path of the file to open

    Returns
    -------
    imps : ImagePlus
        Image opened via BF
    """
    options = ImporterOptions()
    options.setId(str(indivFile))
    options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE)
    imps    = BF.openImagePlus(options)
    return imps


def BFExport(imp, savepath):
    """Export using BioFormats

    Parameters
    ----------
    imp : {ImagePlus}
        ImagePlus of the file to save
    savepath : {str}
        Path where to save the image

    """

    print('Savepath: ', savepath)
    plugin     = LociExporter()
    plugin.arg = savepath
    exporter   = Exporter(plugin, imp)
    exporter.run()

# ─── VARIABLES ──────────────────────────────────────────────────

# Minimum volume to be considered Golgi
volMin  = 5
# Filter objects touching in Z ?
filter_objects_touching_z = False

# OMERO server info
HOST = "xxxx.xxxx.xxxx.xx"
PORT = 4064
datasetId = datasetid
groupId = "-1"

# ─── MAIN CODE ──────────────────────────────────────────────────────────────────

# Prototype analysis example
gateway = omeroConnect()
imageIds = getImageIds(gateway, datasetId)
imageIds.sort()
# reload(sys)
# sys.setdefaultencoding('utf8')

for imageId in imageIds:

    openImagePlus(HOST, USERNAME, PASSWORD, groupId, imageId)
    imp = IJ.getImage()
    IJ.log("Now processing " + imp.getTitle())
    
    origin_folder = str(destination)
    origin_name   = imp.getTitle().replace(" ", "_")
    short_name    = os.path.splitext(origin_name)[0]
    
    out_full_path = origin_folder + os.sep + short_name
    
    channel_of_interest = 2
    
    imp_channel_of_interest = Duplicator().run(imp, channel_of_interest,
                                               channel_of_interest, 1, imp.getNSlices(), 1, 1)
    
    IJ.log("Thresholding...")
    IJ.setAutoThreshold(imp_channel_of_interest, "Otsu dark stack")
    IJ.run(imp_channel_of_interest, "Convert to Mask",
           "method=Otsu background=Dark black")

    imp_channel_of_interest.setCalibration(imp.getCalibration())


    IJ.log("3D segmentation...")
    segment_3D  = Segment3DImage(imp_channel_of_interest, 0, 255)
    segment_3D.segment()
    stack_label = segment_3D.getLabelledObjectsStack()
    imp_label   = ImagePlus("3D Labelled", stack_label)
    imp_label.setCalibration(imp.getCalibration())
    
    # wrap ImagePlus into 3D suite image format
    img  = ImageInt.wrap(imp_label)
    # create a population of 3D objects
    pop  = Objects3DPopulation(img)
    nb   = pop.getNbObjects()
    # print nb
    unit = imp.getCalibration().getUnits()
    # print(unit)
    # print(nb)
    IHimp2              = ImageHandler.wrap(imp)
    volume_list         = []
    # convex_volume_list  = []
    # solidity_3D_list    = []
    mean_intensity_list = []
    feret_list          = []
    compactness_list    = []
    object_number_list  = []
    surface_list        = []
    golgi_obj_to_remove = []

    IJ.log("Looping over the objects...")
    # loop over the objects
    for i in range(0, nb):
        # IJ.log("Object " + str(i+1) + "/" + str(nb))
        obj = pop.getObject(i)
    
        # Check if we should filter out the objects touching in Z
        if(obj.touchBorders(img, filter_objects_touching_z)):
            golgi_obj_to_remove.append(obj)
            continue
        object_volume = obj.getVolumeUnit()
    
        # Discard objects with a volume below 5000µm3
        if(object_volume < volMin):
            golgi_obj_to_remove.append(obj)
            continue
    
        # Get surface
        surface_list.append(obj.getAreaUnit())
        volume_list.append(object_volume)
    
        # Get object number
        object_number_list.append(i)
        
        # # Get the convex hull object
        # IJ.log("Get convex hull object... (takes time)")
        # conv_obj = obj.getConvexObject()
    
    
        # Measure Convex Hull volume
        # convex_volume = conv_obj.getVolumeUnit()
        # convex_volume_list.append(convex_volume)
        # solidity_3D_list.append(object_volume / convex_volume)
    
        # Measure mean intensity
        mean_intensity_list.append(obj.getPixMeanValue(IHimp2))
    
        # Measure feret
        feret_list.append(obj.getFeret())
    
        # Measure compactness
        compactness_list.append(obj.getCompactness(True))

    IJ.log("Saving results...")
    outCSV = out_full_path + ".csv"

    for obj in golgi_obj_to_remove:
        pop.removeObject(obj)
    pop.saveObjects(out_full_path + "_3DROI.zip")

    # imp.close()

    # print(volume_list)
    with open(outCSV, 'wb') as f:
        writer = csv.writer(f)
        # writer.writerow(
        #     ["Object number, Object Volume (" + unit + " cube), Convex Hull Volume (" + unit + " cube), solidity3D, Compactness, Surface (" + unit + " cube), Mean Intensity, Feret diameter (" + unit + ")"])
        writer.writerow(
            ["Object number, Object Volume (" + unit + " cube), Compactness, Surface (" + unit + " cube), Mean Intensity, Feret diameter (" + unit + ")"])
        writer.writerows(
            izip(object_number_list, volume_list, compactness_list, surface_list, mean_intensity_list, feret_list))

    IJ.log("Results for " + origin_name + " have been saved as " + outCSV)
    imp.close()

IJ.log("##############")
IJ.log("MACRO FINISHED")
IJ.log("##############")
