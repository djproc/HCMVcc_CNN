#writing functions to be called from different notebooks, making the code easier to read
import os
import glob
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import shutil
import pandas as pd
import cv2
from tqdm import tqdm
import pickle
import skimage.io
import skimage.measure
import tifffile
import scipy.ndimage
from skimage import exposure

def DataCheck(df_combined_data, group_data_by, display_x_axis, display_y_axis):
    fig, ax =plt.subplots(1,3, figsize=(20,4))

    sns.violinplot(x=group_data_by, y='AreaShape_Area', data=df_combined_data.reset_index(), size=1,ax=ax[0])
    sns.violinplot(x=group_data_by, y=display_y_axis, data=df_combined_data.reset_index(), size=1,ax=ax[1])
    sns.scatterplot(x=display_x_axis, y=display_y_axis, data=df_combined_data.reset_index(), hue=group_data_by, size=0.0001, alpha=0.01,ax=ax[2],legend=False)
    fig.show()
    #TODO: RENDER A FEW CELLS HERE ALSO, RANDOMLY SO WE GET AN IDEA OF DATA QUALITY BEFORE EXPORT

def ImportData_NUC_CYTO():
    cp_output_dir = "_CellProfiler"
    os.makedirs(cp_output_dir, exist_ok=True)
    #print(f'CREATED FOLDER NAMED: {cp_output_dir} \nFOLDER LOCATED AT: {os.getcwd()}')
    
    #TODO: get .csv files and add to this folder, also add any .cppipe file and any file with the same name
    #exp_name = os.path.basename(glob.glob('*_Image.csv',recursive=True)[0])[:-10] #_Image.csv files will always be saved
    
    for filename in Path(os.getcwd()).rglob('*_Image.csv'):
        exp_name = os.path.basename(filename)[:-10]
    pickle.dump( exp_name, open( "HCMVcc_variables.p", "wb" ) )
    print(f'EXPERIMENT NAME: {exp_name}')
    for filename in Path(os.getcwd()).rglob(f'{exp_name}*.csv'):
        #if "linescans" not in filename:
        shutil.move(filename, os.path.join(cp_output_dir,os.path.basename(filename)))
    
    #move files with exp_name in them to the cp_output_dir
    #exp_files = glob.glob(f'{exp_name}*',recursive=True)
    #for f in exp_files:
    #    shutil.move(f, os.path.join(cp_output_dir,f))
        
    #generate the other cellprofiler output filenames
    nuc_csv =  f"{exp_name}_NUC_DAPI.csv"
    cyto_csv = f"{exp_name}_Cytoplasm.csv"
    image_csv = f"{exp_name}_Image.csv"
    print("\nIMPORTED AND MERGED THE FOLLOWING FILES:", nuc_csv, cyto_csv, image_csv, sep="\n - ")
    
    #import these files as datafames using pandas
    #nucleus data
    df_nuc = pd.read_csv(os.path.join(cp_output_dir, nuc_csv), na_filter=True)
    df_nuc.set_index("ImageNumber", inplace=True)
    #cytoplasm data
    df_cyto = pd.read_csv(os.path.join(cp_output_dir, cyto_csv), na_filter=True)
    df_cyto.set_index("ImageNumber", inplace=True)
    #image info
    df_image = pd.read_csv(os.path.join(cp_output_dir, image_csv), na_filter=True)
    df_image.set_index("ImageNumber", inplace=True)
    #then extract only the image urls from this 
    df_image_url = df_image.filter(regex=r'^URL_', axis=1) #this will select any columns starting with "URL_"
    
    #combine these dataframes together
    #merge nucleus data with urls
    df_combined_data = df_nuc.merge(df_image_url, left_on='ImageNumber', right_on='ImageNumber', how='outer')
    #merge this with cytoplasm data and differentiate columns from the two datasets as "_NUC" and "_CYTO"
    df_combined_data = df_combined_data.merge(df_cyto, left_on=["ImageNumber", "ObjectNumber"], right_on=["ImageNumber", "Parent_NUC_DAPI"], how="outer", suffixes=('_NUC', '_CYTO'))

    #we can also just look at the raw number of rows in the dataframe to see how many nuclei we've identified
    df_combined_data.reset_index(inplace=True)
    df_combined_data.rename(columns={'ObjectNumber_NUC':'NUC_ID'}, inplace=True)
    df_combined_data.dropna(subset=['NUC_ID'], inplace=True)
    print(f'\nDETECTED NUCLEI: {df_combined_data.shape[0]:,.0f}')
    return(df_combined_data, df_image_url, exp_name);


def GenerateIDs_IMGexport(df_combined_data, group_data_by,C1,C2,C3,exp_name, export_dimensions):
    
    #conversion to string for concatenation
    df_combined_data["NUC_ID"] = df_combined_data["NUC_ID"].astype(float).astype(int).astype(str)
    df_combined_data["ImageNumber"] = df_combined_data["ImageNumber"].astype(float).astype(int).astype(str)
    df_combined_data['Metadata_date_NUC'] = df_combined_data['Metadata_date_NUC'].astype(float).astype(int).astype(str)
    df_combined_data['Metadata_experiment_NUC'] = df_combined_data['Metadata_experiment_NUC'].astype(float).astype(int).astype(str) 
    # TODO: if this is a string and not a number, then make it an float then int then string. 
    df_combined_data['Metadata_for_ID'] = df_combined_data[group_data_by].astype(str)
    # If there is no biorep value, set it to 'repX'
    if 'Metadata_biorep_NUC' not in df_combined_data:
        df_combined_data['Metadata_biorep_NUC'] = 'repX' 
    df_combined_data["Unique_ID"] = df_combined_data['Metadata_experiment_NUC'].str.cat(df_combined_data[['Metadata_date_NUC', 
                                                                                                          'Metadata_biorep_NUC',
                                                                                                          'Metadata_for_ID',
                                                                                                          'ImageNumber','NUC_ID']
                                                                                                        ], sep="_")
    print(f'EXAMPLE IDS: {df_combined_data["Unique_ID"][1]}')
    
    img_C1_16bit = skimage.io.imread(df_combined_data[C1][1], as_gray=True, plugin='tifffile')
    original_img_dimensions = img_C1_16bit.shape[0]
    cp_downscaled_ratio = original_img_dimensions / 256 # CellProfiler pipeline downscales images to 256x256 for increased speed
    print(f'CellProfiler scaled images by {cp_downscaled_ratio}')
    img_C1_16bit = None
    
    df_combined_data['NUC_x0'] = df_combined_data[('Location_Center_X_NUC')] * cp_downscaled_ratio #correcting for downscaling
    df_combined_data['NUC_y0'] = df_combined_data[('Location_Center_Y_NUC')] * cp_downscaled_ratio #correcting for downscaling
    
    df_combined_data['CYTO_x0'] = df_combined_data[('Location_CenterMassIntensity_X_gB_small')] * cp_downscaled_ratio #correcting for downscaling
    df_combined_data['CYTO_y0'] = df_combined_data[('Location_CenterMassIntensity_Y_gB_small')] * cp_downscaled_ratio #correcting for downscaling
    
    df_combined_data['URL_C1'] = df_combined_data[C1] #this will be red in imageJ
    df_combined_data['URL_C2'] = df_combined_data[C2] #this will be green in imageJ
    df_combined_data['URL_C3'] = df_combined_data[C3] #this will be blue in imageJ
    
    df_display = df_combined_data[["Unique_ID",
                                   'NUC_x0',
                                   'NUC_y0',
                                   'URL_C1',
                                   'URL_C2',
                                   'URL_C3',
                                   'AreaShape_Orientation',
                                   'CYTO_x0',
                                   'CYTO_y0']]
        
    #create a folder to export data from the CNN
    cnn_export_folder = f"_IMGexportforCNN_{exp_name}"
    os.makedirs(cnn_export_folder, exist_ok=True)

    #just in case this crashes half way through, if you re-run this function this checks what you have already exported and continues from there
    already_exported_files = [os.path.splitext(filename)[0] for filename in os.listdir(cnn_export_folder)]
    
    #then run through the dataframe for all values that haven't been exported yet
    for index, row in tqdm(df_display.iterrows(), total=df_display.shape[0]):
        if df_display['Unique_ID'][index] not in already_exported_files:
    
            export_dimensions = int(export_dimensions)
            img_border = int(export_dimensions/2)
            dpi = int(export_dimensions/2)
            
            img_C1_16bit = skimage.io.imread(df_combined_data["URL_C1"][index], as_gray=True, plugin='tifffile')
            img_C2_16bit = skimage.io.imread(df_combined_data["URL_C2"][index], as_gray=True, plugin='tifffile')
            #img_C3_16bit = skimage.io.imread(df_combined_data["URL_C3"][index], as_gray=True, plugin='tifffile')
            #img_C4_16bit = skimage.io.imread(df_combined_data["URL_C4"][index], as_gray=True, plugin='tifffile')

            p0, p100 = np.percentile(img_C1_16bit, (0, 100))
            img_C1_16bit_rescaled = exposure.rescale_intensity(img_C1_16bit, in_range=(p0*1.3, p100*1.1))

            p0, p100 = np.percentile(img_C2_16bit, (0, 100))
            img_C2_16bit_rescaled = exposure.rescale_intensity(img_C2_16bit, in_range=(p0*1.3, p100))

            img_C1_8bit = (img_C1_16bit_rescaled/256).astype('uint8')
            img_C2_8bit = (img_C2_16bit_rescaled/256).astype('uint8')
            img_blank_8bit = img_C1_8bit * 0 # this makes it black
            #img_C3_8bit = img_C3_16bit/32).astype('uint8')
            #img_C4_8bit = img_C3_16bit/32).astype('uint8')

            #to save on memory, now that we've used this data, set it to None

            merged_channels = cv2.merge((img_C1_8bit,img_C2_8bit,img_blank_8bit))
            merged_channels_border = cv2.copyMakeBorder(merged_channels, img_border, img_border, img_border, img_border, cv2.BORDER_CONSTANT, value=0)

            x0 = df_display['NUC_x0'][index]
            y0 = df_display['NUC_y0'][index]

            fig = plt.figure()
            fig.set_size_inches(2,2)
            ax = plt.Axes(fig, [0., 0., 1., 1.])
            ax.set_axis_off()
            fig.add_axes(ax)
            plt.xlim(x0,x0+export_dimensions)
            plt.ylim(y0,y0+export_dimensions)
            ax.imshow(merged_channels_border)
            filename = os.path.join(f'{cnn_export_folder}',(df_display['Unique_ID'][index] + '.jpg'))
            fig.savefig(filename, dpi=dpi, bbox_inches='tight', pad_inches=0)

            fig.clf()
            plt.close(fig)
            plt.close('all')
            
            img_C1_16bit = None
            img_C2_16bit = None
            #img_C3_16bit = None
            #img_C4_16bit = None
            img_C1_8bit = None
            img_C2_8bit = None
            img_blank_8bit = None
            merged_channels = None
            merged_channels_border = None

def linescan_calculations(df_combined_data, df_predictions, linescan_len):
    #merge dataframes 
    df_predictions_coords = df_combined_data.merge(df_predictions, left_on='Unique_ID', right_on='Unique_ID', how='right')
    
    #perfrom calculations
    radius_len = 250
    #get the orientation of the minor axis in radians
    
    df_predictions_coords['NUC_orientation_radians'] = (df_predictions_coords[('AreaShape_Orientation')]+90).apply(np.radians)
    df_predictions_coords['NUC_orientation_degrees'] = df_predictions_coords['NUC_orientation_radians'].apply(np.degrees)
    
    #calculate the coordinates of points (x1,y1), (x2,y2) drawing a 1000 pixel line along the minor axis and passing through the center of the nucleus 
    df_predictions_coords['NUC_x1'] = df_predictions_coords['NUC_x0'] + df_predictions_coords['NUC_orientation_radians'].apply(np.cos).multiply(radius_len)
    df_predictions_coords['NUC_y1'] = df_predictions_coords['NUC_y0'] + df_predictions_coords['NUC_orientation_radians'].apply(np.sin).multiply(radius_len)
    
    df_predictions_coords['NUC_x2'] = df_predictions_coords['NUC_x0'] - df_predictions_coords['NUC_orientation_radians'].apply(np.cos).multiply(radius_len)
    df_predictions_coords['NUC_y2'] = df_predictions_coords['NUC_y0'] - df_predictions_coords['NUC_orientation_radians'].apply(np.sin).multiply(radius_len)
    
    
    #THE FOLLOWING CALCULATIONS ARE ONLY RELEVANT FOR TB96 BUT WILL ALSO BE PERFORMED ON MOCK
    #get the center of the AC (AC_x0)(we will filter out for uninfected cells later)
    df_predictions_coords['AC_x0'] = df_predictions_coords['CYTO_x0']
    df_predictions_coords['AC_y0'] = df_predictions_coords['CYTO_y0']
    
    #calculate the AC to NUC orientation in radians 
    df_predictions_coords["AC_to_NUC_orientation_radians"] =  np.arctan2((df_predictions_coords['AC_y0'] - df_predictions_coords['NUC_y0']),(df_predictions_coords['AC_x0'] - df_predictions_coords['NUC_x0']))
    df_predictions_coords.loc[df_predictions_coords["AC_to_NUC_orientation_radians"]<0, ["AC_to_NUC_orientation_radians"]] = 2*np.pi+ df_predictions_coords["AC_to_NUC_orientation_radians"]
    
    df_predictions_coords["ACtoNUC_X_Midpoint"] = (df_predictions_coords['NUC_x0'] + df_predictions_coords['AC_x0'])/2
    df_predictions_coords["ACtoNUC_Y_Midpoint"] = (df_predictions_coords['NUC_y0'] + df_predictions_coords['AC_y0'])/2
    
    df_predictions_coords["X_AC_fixedradius"] = df_predictions_coords["ACtoNUC_X_Midpoint"] + radius_len *  np.cos(df_predictions_coords['AC_to_NUC_orientation_radians'] - np.pi)
    df_predictions_coords["Y_AC_fixedradius"] = df_predictions_coords["ACtoNUC_Y_Midpoint"] + radius_len *  np.sin(df_predictions_coords['AC_to_NUC_orientation_radians'] - np.pi)
    
    df_predictions_coords["X_NUC_fixedradius"] = df_predictions_coords["ACtoNUC_X_Midpoint"] + radius_len * np.cos(df_predictions_coords['AC_to_NUC_orientation_radians'])
    df_predictions_coords["Y_NUC_fixedradius"] = df_predictions_coords["ACtoNUC_Y_Midpoint"] + radius_len * np.sin(df_predictions_coords['AC_to_NUC_orientation_radians'])
    
    df_predictions_coords["AC_to_NUC_orientation_degrees"] = df_predictions_coords['AC_to_NUC_orientation_radians'].apply(np.degrees)
    
    return(df_predictions_coords)

class data_linewidth_plot():
    def __init__(self, x, y, **kwargs):
        self.ax = kwargs.pop("ax", plt.gca())
        self.fig = self.ax.get_figure()
        self.lw_data = kwargs.pop("linewidth", 1)
        self.lw = 1
        self.fig.canvas.draw()

        self.ppd = 72./self.fig.dpi
        self.trans = self.ax.transData.transform
        self.linehandle, = self.ax.plot([],[],**kwargs)
        if "label" in kwargs: kwargs.pop("label")
        self.line, = self.ax.plot(x, y, **kwargs)
        self.line.set_color(self.linehandle.get_color())
        self._resize()
        self.cid = self.fig.canvas.mpl_connect('draw_event', self._resize)

    def _resize(self, event=None):
        lw =  ((self.trans((1, self.lw_data))-self.trans((0, 0)))*self.ppd)[1]
        if lw != self.lw:
            self.line.set_linewidth(lw)
            self.lw = lw
            self._redraw_later()

    def _redraw_later(self):
        self.timer = self.fig.canvas.new_timer(interval=10)
        self.timer.single_shot = True
        self.timer.add_callback(lambda : self.fig.canvas.draw_idle())
        self.timer.start()

def line_profile_coordinates(src, dst, linewidth=1):
    """
    TAKEN FROM: https://github.com/scikit-image/scikit-image/blob/master/skimage/measure/profile.py
    Return the coordinates of the profile of an image along a scan line.
    Parameters
    ----------
    src : 2-tuple of numeric scalar (float or int)
        The start point of the scan line.
    dst : 2-tuple of numeric scalar (float or int)
        The end point of the scan line.
    linewidth : int, optional
        Width of the scan, perpendicular to the line
    Returns
    -------
    coords : array, shape (2, N, C), float
        The coordinates of the profile along the scan line. The length of the
        profile is the ceil of the computed length of the scan line.
    Notes
    -----
    This is a utility method meant to be used internally by skimage functions.
    The destination point is included in the profile, in contrast to
    standard numpy indexing.
    """
    src_row, src_col = src = np.asarray(src, dtype=float)
    dst_row, dst_col = dst = np.asarray(dst, dtype=float)
    d_row, d_col = dst - src
    theta = np.arctan2(d_row, d_col)

    length = int(np.ceil(np.hypot(d_row, d_col) + 1))
    # we add one above because we include the last point in the profile
    # (in contrast to standard numpy indexing)
    line_col = np.linspace(src_col, dst_col, length)
    line_row = np.linspace(src_row, dst_row, length)

    # we subtract 1 from linewidth to change from pixel-counting
    # (make this line 3 pixels wide) to point distances (the
    # distance between pixel centers)
    col_width = (linewidth - 1) * np.sin(-theta) / 2
    row_width = (linewidth - 1) * np.cos(theta) / 2
    perp_rows = np.array([np.linspace(row_i - row_width, row_i + row_width,
                                      linewidth) for row_i in line_row])
    perp_cols = np.array([np.linspace(col_i - col_width, col_i + col_width,
                                      linewidth) for col_i in line_col])
    return np.array([perp_rows, perp_cols])

def getlinescans_MOCK(df_predictions_coords_MOCK, channel_li, exp_name, linewidth):
    
    #create a new dataframe to add the linescan values to 
    df_linescans_MOCK = pd.DataFrame()

    #then loop through the indexes in the dataframe with the coordinates (note tqdm just adds a timer to this function)
    for index, row in tqdm(df_predictions_coords_MOCK.iterrows(), total=df_predictions_coords_MOCK.shape[0]):
        
        x1 = df_predictions_coords_MOCK["NUC_x1"][index]
        y1= df_predictions_coords_MOCK["NUC_y1"][index]
    
        x2 = df_predictions_coords_MOCK["NUC_x2"][index]
        y2 = df_predictions_coords_MOCK["NUC_y2"][index]
        
        for c in channel_li:
            ID_Channel = df_predictions_coords_MOCK["Unique_ID"][index] + "_" + c
            
            img_url = df_predictions_coords_MOCK[c][index]
            #image = skimage.io.imread(df_predictions_coords[c][index], as_gray=True)
            img = skimage.io.imread(img_url, as_gray=True, plugin='tifffile')
            #linescan = pd.Series(skimage.measure.profile_line(image, src, dst, linewidth=150, order=0, mode='constant', cval=0.0, reduce_func=np.mean))
            # i used this over scikit-image as there was a bug in their function, this should be more stable
            linescan = scipy.ndimage.map_coordinates(np.transpose(img), line_profile_coordinates((x1,y1),(x2,y2), linewidth=linewidth))
            linescan = np.ma.masked_equal(linescan,0)
            linescan_mean = linescan.mean(axis=1)
    
            df_linescans_MOCK[ID_Channel] = pd.Series(linescan_mean)
    
    #df_linescans_MOCK.replace(0, np.nan, inplace=True)
    os.makedirs("_linescans", exist_ok=True)
    df_linescans_MOCK.to_csv(f'_linescans/{exp_name}_linescans_MOCK.csv',index=False)
    return(df_linescans_MOCK)

def getlinescans_TB96(df_predictions_coords_TB96, channel_li, exp_name, linewidth):
    
    #create a new dataframe to add the linescan values to 
    df_linescans_TB96 = pd.DataFrame()

    #then loop through the indexes in the dataframe with the coordinates (note tqdm just adds a timer to this function)
    for index, row in tqdm(df_predictions_coords_TB96.iterrows(), total=df_predictions_coords_TB96.shape[0]):
         
        x1 = df_predictions_coords_TB96["X_AC_fixedradius"][index]
        y1 = df_predictions_coords_TB96["Y_AC_fixedradius"][index]
        
        x2 = df_predictions_coords_TB96["X_NUC_fixedradius"][index]
        y2 = df_predictions_coords_TB96["X_NUC_fixedradius"][index]
        
        for c in channel_li:
            ID_Channel = df_predictions_coords_TB96["Unique_ID"][index] + "_" + c
            
            img_url = df_predictions_coords_TB96[c][index]
            #image = skimage.io.imread(df_predictions_coords[c][index], as_gray=True)
            img = skimage.io.imread(img_url, as_gray=True, plugin='tifffile')
            #linescan = pd.Series(skimage.measure.profile_line(image, src, dst, linewidth=150, order=0, mode='constant', cval=0.0, reduce_func=np.mean))
            # i used this over scikit-image as there was a bug in their function, this should be more stable
            linescan = scipy.ndimage.map_coordinates(np.transpose(img), line_profile_coordinates((x1,y1),(x2,y2), linewidth=linewidth))
            linescan = np.ma.masked_equal(linescan,0)
            linescan_mean = linescan.mean(axis=1)
    
            df_linescans_TB96[ID_Channel] = pd.Series(linescan_mean)
    
    #df_linescans_MOCK.replace(0, np.nan, inplace=True)
    os.makedirs("_linescans", exist_ok=True)
    df_linescans_TB96.to_csv(f'_linescans/{exp_name}_linescans_TB96.csv',index=False)
    return(df_linescans_TB96)

def preview_linescans_MOCK(df_predictions_coords, img_val, channel_li, linewidth):
    
    channel_li.append("MERGE")
    df_predictions_coords_MOCK = df_predictions_coords[df_predictions_coords["Group"]=="MOCK"]
    
    #TODO: turn this into a for loop over x random examples

    #get coordinates caluclated for (x1,y1) and (x2,y2) to draw linescans
    x1 = df_predictions_coords_MOCK["NUC_x1"][img_val]
    y1 = df_predictions_coords_MOCK["NUC_y1"][img_val]
    x2 = df_predictions_coords_MOCK["NUC_x2"][img_val]
    y2 = df_predictions_coords_MOCK["NUC_y2"][img_val]

    #load images
    img_r = skimage.io.imread(df_predictions_coords_MOCK["URL_C1"][img_val], as_gray=True, plugin='tifffile')
    img_g = skimage.io.imread(df_predictions_coords_MOCK["URL_C2"][img_val], as_gray=True, plugin='tifffile')
    img_b = skimage.io.imread(df_predictions_coords_MOCK["URL_C3"][img_val], as_gray=True, plugin='tifffile')

    #merge rgb
    img_rgb = (np.dstack((img_r, img_g, img_b))/(32)).astype(np.uint8)
    
    #collect linescan data for lines across the full specified line width
    linescan_r = scipy.ndimage.map_coordinates(np.transpose(img_r), line_profile_coordinates((x1,y1),(x2,y2), linewidth=linewidth))
    linescan_b = scipy.ndimage.map_coordinates(np.transpose(img_g), line_profile_coordinates((x1,y1),(x2,y2), linewidth=linewidth))
    linescan_g = scipy.ndimage.map_coordinates(np.transpose(img_b), line_profile_coordinates((x1,y1),(x2,y2), linewidth=linewidth))
    
    #mask zero values within linescan data
    linescan_r = np.ma.masked_equal(linescan_r,0)
    linescan_g = np.ma.masked_equal(linescan_g,0)
    linescan_b = np.ma.masked_equal(linescan_b,0)
    
    #collect the mean linescan value across all lines in width
    linescan_r_mean = linescan_r.mean(axis=1)
    linescan_g_mean = linescan_g.mean(axis=1)
    linescan_b_mean = linescan_b.mean(axis=1)
    
    fig, ax = plt.subplots(2,4, figsize=(20,10), constrained_layout=True)
    
    #display individual channel images and merged image
    ax[0,0].imshow(img_r, cmap='jet')
    ax[0,1].imshow(img_g, cmap='jet')
    ax[0,2].imshow(img_b, cmap='jet')
    ax[0,3].imshow(img_rgb)
    
    #display linescans individually
    ax[1,0].plot(linescan_r_mean)
    ax[1,1].plot(linescan_g_mean)
    ax[1,2].plot(linescan_b_mean)
    
    #display linescans together below the merged iamge
    ax[1,3].plot(linescan_r_mean)
    ax[1,3].plot(linescan_g_mean)
    ax[1,3].plot(linescan_b_mean)                                         
    
    #add titles to each channel column
    for i in range(len(channel_li)):
        ax[0,i].set_title(f'{channel_li[i]}')
    
    #collect other data for drawing annotations on graphs
    midpoint = 500
    span = 100
    img_size = 600
    x_mid = df_predictions_coords_MOCK["NUC_x0"][img_val]
    y_mid = df_predictions_coords_MOCK["NUC_y0"][img_val]
    
    for i in [0,1,2,3]:

        #draw the linescan on the images
        x_scatter = (x_mid,x1,x2)
        y_scatter = (y_mid,y1,y2) 
        l = data_linewidth_plot(x_scatter, y_scatter, ax=ax[0,i], linewidth=linewidth, alpha=0.25, color='white')
        ax[0,i].scatter(x_scatter, y_scatter, s=100,color=('blue','green','red'))
        labels=['mid','source','dest']
        for label_num, label_text in enumerate(labels):
            ax[0,i].annotate(label_text, (x_scatter[label_num], y_scatter[label_num]))

        #draw verticle lines on the graph
        ax[1,i].axvline(x=midpoint, dashes=[6,6], linewidth=1,color='red',alpha=0.5)
        ax[1,i].axvline(x=midpoint-span, dashes=[6,6], linewidth=1,color='grey',alpha=0.5)
        ax[1,i].axvline(x=midpoint+span, dashes=[6,6], linewidth=1,color='grey',alpha=0.5)
        ax[1,i].axvline(x=midpoint-2*span, dashes=[6,6], linewidth=1,color='grey',alpha=0.5)
        ax[1,i].axvline(x=midpoint+2*span, dashes=[6,6], linewidth=1,color='grey',alpha=0.5)
        ax[1,i].axvspan((midpoint-span), midpoint, alpha=0.05, color='grey')
        ax[1,i].axvspan(midpoint,(midpoint+span), alpha=0.05, color='grey')
        
        #set x and y limits for the images 
        ax[0,i].set_xlim(x_mid-img_size,x_mid+img_size)
        ax[0,i].set_ylim(y_mid+img_size,y_mid-img_size) # need to set it like this or otherwise upside down
        
        #set x limis for the linescan graphs
        ax[1,i].set_xlim(midpoint - img_size/2, midpoint + img_size/2)
        
        #remove ticks and labels from some figures
        ax[1,i].set_yticklabels([])
        ax[1,i].set_yticks([])
        ax[0,i].set_yticklabels([])
        ax[0,i].set_yticks([])
        ax[0,i].set_xticklabels([])
        ax[0,i].set_xticks([])

    ax[1,3].yaxis.tick_right()

    plt.title(df_predictions_coords_MOCK["Unique_ID"][img_val])
    plt.show()

def merge_MOCK_TB(df_linescans_MOCK, df_linescans_TB96):
    if df_linescans_MOCK is None:
        df_linescans_MOCK = pd.read_csv(f'_linescans/{exp_name}_linescans_MOCK.csv')
    if df_linescans_TB96 is None:
        df_linescans_TB96 = pd.read_csv(f'_linescans/{exp_name}_linescans_TB96.csv') 
    df_FINAL = pd.concat([df_linescans_MOCK, df_linescans_TB96], axis=1)
    return(df_FINAL)

    
#### BELOW ARE FUNCTIONS YOU CAN RUN IF YOU ANNOTATED YOUR DATA INCORRECTLY, IT HAPPENS... ###
    
def URLFIX(df, find_string, replace_string, columns):
    """Use this to fix path URL issues between operating systems"""
    for col in columns:
        df[col] = df[col].str.replace(find_string,replace_string)
        print(df[col][1])
    return(df)

def colSWAP(df, col_name_A, col_name_B):
    """If you've named your metadata incrorrectly in CellProfiler, you can easily swap dataframe column names with this function"""
    df = df.rename(columns={col_name_A:col_name_B, col_name_B:col_name_A})
    return(df)




# TODO: linescan to analysis functions


"""
def LoadIDsfromCNN(exp_name, date_var, time_var, folder_name_1, folder_name_2):
    #exp_var = "137_20190805_SUN2quant_rerpCGH_v1"
    #date_var = 20190817
    #time_var = "1145"
    #### make this a drop down box ####
    #MOCK_confidece_list = [99,90]
    #TB_confidence_list = [99,90]

    FOLDER_1_files = []
    for conf in FOLDER_1_confidece_list:
        FOLDER_1_files_single_folder = os.listdir(f'{PATH}\\CNNpredictions_{date_var}_{time_var}\\{conf}confidence\\{folder_name_1}')
        print(f'Loading {len(MOCK_files_single_folder)} files for MOCK {conf}confidence')
        MOCK_files = MOCK_files + MOCK_files_single_folder
    print(f'Loaded {len(MOCK_files)} files total for TB')

    TB_files = []
    for conf in TB_confidence_list:
        TB_files_single_folder = os.listdir(f'{PATH}\\CNNpredictions_{date_var}_{time_var}\\{conf}confidence\\TB_perfect')
        print(f'Loading {len(TB_files_single_folder)} files for TB {conf}confidence')
        TB_files = TB_files + TB_files_single_folder
    print(f'Loaded {len(TB_files)} files total for TB')

    #convert these lists to a dataframe
    df_TB_ID = pd.DataFrame({'Unique_ID':TB_files})
    df_MOCK_ID = pd.DataFrame({'Unique_ID':MOCK_files})

    #strip "_RGB.jpg" from filename to match UNIQUE ID
    df_TB_ID['Unique_ID'].replace(regex=True,inplace=True,to_replace='_RGB.jpg',value='')
    df_MOCK_ID['Unique_ID'].replace(regex=True,inplace=True,to_replace='_RGB.jpg',value='')

    #filter mock for only in mock sample (there are obviously some uninfected cells in the TB infected samples)
    df_MOCK_ID = df_MOCK_ID[df_MOCK_ID['Unique_ID'].str.contains('MOCK')]
    print(df_MOCK_ID.shape)
    #filter TB for only in TB sample (just in case)
    df_TB_ID = df_TB_ID[df_TB_ID['Unique_ID'].str.contains('_96hpi_')]
    print(df_TB_ID.shape)

    # IMPORT DATAFRAME WHERE WE CREATED THE IDs and COORDS FOR THIS DATASET
    df_coords = pd.read_csv(f'{exp_var}_IDsandCoords.csv')
    #df_FULL = pd.read_csv("117_20181023_SUN2_v4_IDsandCoords_FULL.csv") #create a FULL version of the dataframe in the Generate IDs and coords notebook

    # Merge coords with lists from CNN sorting for TB
    df_TB_CNNsorted99_coords = df_coords.merge(df_TB_ID, left_on='Unique_ID', right_on='Unique_ID', how='right')
    df_MOCK_CNNsorted99_coords = df_coords.merge(df_MOCK_ID, left_on='Unique_ID', right_on='Unique_ID', how='right')
"""

# TODO: merging and graphing functions


