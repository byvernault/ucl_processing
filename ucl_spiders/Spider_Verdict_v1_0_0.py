"""Spider_Verdict.

Author:         Benjamin Yvernault
contact:        b.yvernault@ucl.ac.uk
Spider name:    Verdict
Spider version: 1.0.0
Creation date:  2016-08-10 17:48:51.701453
Purpose:        Generate Verdict Map from all Verdict scans registered \
                together and merge as one nifti
"""

# Python packages import
import os
import sys
import dicom
import shutil
import datetime
import numpy as np
import nibabel as nib
from dicom.sequence import Sequence
from dicom.dataset import Dataset, FileDataset
from dax import XnatUtils, spiders, SessionSpider

__author__ = "Benjamin Yvernault"
__email__ = "b.yvernault@ucl.ac.uk"
__purpose__ = "Generate Verdict Map from all Verdict scans registered together \
and merge as one nifti"
__spider_name__ = "Verdict"
__version__ = "1.0.0"
__modifications__ = """2016-08-10 17:48:51.701453 - Original write
2016-12-13 11:21:30 - Re-organisation of the matlab folder. Removing amico \
options."""


DEFAULT_MODEL = "VerdictProstate_Rmaps"
DEFAULT_SCHEME_FILE = "NOptimisedV_IN.scheme"
DEFAULT_VERDICT_TEMPLATE = """
addpath(genpath('{matlab_code}'));
launch_AMICO_for_INNOVATE('{input_path}',\
'{subject}',\
'{filename}',\
'{output}',\
'{project}',\
'{matlab_code}/AMICO/matlab/',\
'{camino}',\
'{spams}',\
'{scheme_filename}',\
'{model}');
"""
DEFAULT_PDF_MAKER = """
% Code to generate PDF for VERDICT MAP

% Add path where matlab Code is
addpath(genpath('{matlab_code}'));

% Maps Name:
nb_acq = {acq};
maps_name = {{'fIC','R','cellularity','fEES','fVASC','FobjCamino'}};

% Open an image to see the number of slices
f_file = fullfile('{maps_folder}', 'FIT_FobjCamino.nii');
aux = load_untouch_nii(f_file);
nb_slices = aux.hdr.dime.dim(4);
for i=1:nb_slices
    filename = ['DisplayMapsVerdictSlice' num2str(i, '%02d');];
    jpg_path = fullfile('{output_folder}', [filename '.pdf']);
    plot_oneslice_selectedmaps('{maps_folder}','{subject}',i,maps_name,1,1,1,jpg_path,nb_acq);
end
"""
DICOM_SCAN_TYPE = ['WIP b3000_90 SENSE', 'SWITCH DB TO YES b3000_80',
                   'b3000_80']
ORDER_MAPS = ['FIT_fIC', 'FIT_cellularity',
              'FIT_fEES', 'FIT_fVASC', 'FIT_R',
              'FIT_FobjCamino', 'FIT_Fobj', 'FIT_dir',
              'FIT_R_0-3u', 'FIT_R_3-6u', 'FIT_R_6-9u',
              'FIT_R_9-12u', 'FIT_R_12-15u',
              'FIT_fIC_0-3u', 'FIT_fIC_3-6u', 'FIT_fIC_6-9u',
              'FIT_fIC_9-12u', 'FIT_fIC_12-15u',
              'FIT_cell_0-3u', 'FIT_cell_3-6u', 'FIT_cell_6-9u',
              'FIT_cell_9-12u', 'FIT_cell_12-15u']
C_RANGE = {'FIT_fIC': {'min': 0, 'max': 1, 'index': 1},
           'FIT_cellularity': {'min': 2*10**11, 'max': 1.5*10**14,
                               'index': 3},
           'FIT_fEES': {'min': 0, 'max': 1, 'index': 5},
           'FIT_fVASC': {'min': 0, 'max': 1, 'index': 7},
           'FIT_R': {'min': 0, 'max': 15.10*10**-6, 'index': 9},
           'FIT_FobjCamino': {'min': 0, 'max': 50, 'index': 11},
           'FIT_Fobj': {'min': 0, 'max': 1, 'index': 13},
           'FIT_dir': {'min': 1*10**-10, 'max': 2.9*10**-9, 'index': 15},
           'FIT_R_0-3u': {'min': 0, 'max': 2.67*10**-6, 'index': 17},
           'FIT_R_3-6u': {'min': 3.56*10**-6, 'max': 5.34*10**-6,
                          'index': 19},
           'FIT_R_6-9u': {'min': 6.22*10**-6, 'max': 8.89*10**-6,
                          'index': 21},
           'FIT_R_9-12u': {'min': 9.77*10**-6, 'max': 11.55*10**-6,
                           'index': 23},
           'FIT_R_12-15u': {'min': 12.44*10**-6, 'max': 15.10*10**-6,
                            'index': 25},
           'FIT_fIC_0-3u': {'min': 0, 'max': 1, 'index': 27},
           'FIT_fIC_3-6u': {'min': 0, 'max': 1, 'index': 29},
           'FIT_fIC_6-9u': {'min': 0, 'max': 1, 'index': 31},
           'FIT_fIC_9-12u': {'min': 0, 'max': 1, 'index': 33},
           'FIT_fIC_12-15u': {'min': 0, 'max': 1, 'index': 35},
           'FIT_cell_0-3u': {'min': 3*10**12, 'max': 5*10**16, 'index': 37},
           'FIT_cell_3-6u': {'min': 2*10**12, 'max': 3*10**14, 'index': 39},
           'FIT_cell_6-9u': {'min': 1*10**12, 'max': 2.5*10**14,
                             'index': 41},
           'FIT_cell_9-12u': {'min': 1*10**12, 'max': 5*10**13,
                              'index': 43},
           'FIT_cell_12-15u': {'min': 2*10**11, 'max': 5.5*10**13,
                               'index': 45}
           }


def parse_args():
    """Argument parser for the spider input variables.

    by default (set in get_session_argparser):
        -p       : proj_label
        -s       : subj_label
        -e       : sess_label
        -d       : temp_dir
        --suffix : suffix (suffix for assessor proctype on XNAT)
        --host : host for XNAT (default: XNAT_HOST env variable)
        --user : user on XNAT (default: XNAT_USER env variable)
    your arguments:
        --proctype: proctype where inputs are
        --mc: matlab code to launch verdict map
        --camino : path to Camino
        --spams : path to spams-matlab

    :return: argument parser object created by parse_args()
    """
    ap = spiders.get_session_argparser("Verdict", __purpose__)
    ap.add_argument("--proctype", dest="proctype", default=None, required=True,
                    help="Assessor type containing the registered nifti.")
    ap.add_argument("--nbAcq", dest="nb_acquisition", default=1,
                    required=True,
                    help="Number of Acquisition of VERDICT scans (1 or 2).")
    ap.add_argument("--mc", dest="matlab_code", default=None, required=True,
                    help="Matlab code folder where is \
launch_AMICO_for_INNOVATE and AMICO code in a folder AMICO/matlab.")
    ap.add_argument("--camino", dest="camino", default=None, required=True,
                    help="Path to Camino folder.")
    ap.add_argument("--spams", dest="spams", default=None, required=True,
                    help="Path to spams-matlab folder.")
    ap.add_argument("--scheme", dest="scheme_filename", default=None,
                    help="Path to schemeFilename \
(NOptimisedV_IN.scheme).")
    ap.add_argument("--model", dest="model", default=DEFAULT_MODEL,
                    help="Model type for AMICO: VerdictProstate_Rmaps.")
    return ap.parse_args()


class Spider_Verdict(SessionSpider):
    """Session Spider: Spider_Verdict.

    :param spider_path: spider file path
    :param jobdir: directory for temporary files
    :param xnat_project: project ID on XNAT
    :param xnat_subject: subject label on XNAT
    :param xnat_session: experiment label on XNAT
    :param matlab_path: path to matlab code
    :param xnat_host: host for XNAT if not set in environment variables
    :param xnat_user: user for XNAT if not set in environment variables
    :param xnat_pass: password for XNAT if not set in environment variables
    :param suffix: suffix to the assessor creation
    """

    def __init__(self, spider_path, jobdir, xnat_project, xnat_subject,
                 xnat_session, proctype, nb_acquisition, matlab_code, camino,
                 spams, scheme_filename=None, model=DEFAULT_MODEL,
                 xnat_host=None, xnat_user=None, xnat_pass=None, suffix=""):
        """Entry point for Spider_Verdict Class."""
        super(Spider_Verdict,
              self).__init__(spider_path, jobdir,
                             xnat_project, xnat_subject, xnat_session,
                             xnat_host, xnat_user, xnat_pass,
                             suffix)
        self.proctype = proctype
        self.nb_acquisition = int(nb_acquisition)

        self.inputs = dict()
        self.pdf_final = os.path.join(self.jobdir,
                                      '%s_VERDICT_report.pdf' % xnat_session)

        self.matlab_code = matlab_code
        self.spams = spams
        self.camino = camino
        if not scheme_filename:
            self.scheme_filename = os.path.join(matlab_code,
                                                DEFAULT_SCHEME_FILE)
        else:
            self.scheme_filename = scheme_filename
        self.model = model

    def pre_run(self):
        """Method to download data from XNAT.

        :param argument_parse: argument parser object return by parse_args()
        """
        resource = 'ACQ'  # resource to download from the scan on XNAT
        folder = os.path.join(self.jobdir, 'inputs')
        os.makedirs(folder)
        assessor_label = '-x-'.join([self.xnat_project,
                                     self.xnat_subject,
                                     self.xnat_session,
                                     self.proctype])

        xnat = XnatUtils.get_interface()
        a = xnat.select('/projects/%s/subjects/%s/experiments/%s/assessors/%s'
                        % (self.xnat_project, self.xnat_subject,
                           self.xnat_session, assessor_label))

        for nb_acq in range(1, self.nb_acquisition+1):
            res_name = '%s%d' % (resource, nb_acq)
            self.time_writer('Download resource: %s' % res_name)
            file_gzip = XnatUtils.download_file_from_obj(
                                  folder, a.out_resource(res_name))
            self.time_writer('Unzip file: %s' % file_gzip)
            XnatUtils.gunzip_file(file_gzip)
            self.inputs[nb_acq] = file_gzip[:-3]

        sc = XnatUtils.get_good_scans(a.parent(), DICOM_SCAN_TYPE)[0]
        dcm_file = XnatUtils.download_file_from_obj(
                              folder, sc.resource('DICOM'))
        self.inputs['dcm'] = dcm_file
        xnat.disconnect()

    def run(self):
        """Method running the process for the spider on the inputs data."""
        output_folder = XnatUtils.makedir(os.path.join(self.jobdir, 'outputs'))
        osirix_folder = XnatUtils.makedir(os.path.join(output_folder,
                                                       'OsiriX'))
        for nb_acq in range(1, self.nb_acquisition+1):
            folder = os.path.join(output_folder, str(nb_acq))
            os.makedirs(folder)
            mat_lines = DEFAULT_VERDICT_TEMPLATE.format(
                    matlab_code=self.matlab_code,
                    input_path=os.path.dirname(self.inputs[nb_acq]),
                    subject=self.xnat_subject,
                    filename=os.path.basename(self.inputs[nb_acq]),
                    output=folder,
                    project=self.xnat_project,
                    camino=self.camino,
                    spams=self.spams,
                    scheme_filename=self.scheme_filename,
                    model=self.model)
            matlab_script = os.path.join(output_folder,
                                         'run_verdict_map%d.m' % nb_acq)
            with open(matlab_script, "w") as f:
                f.writelines(mat_lines)
            self.run_matlab(matlab_script, verbose=True)

            # Generate Dicom for OsiriX
            outdir = os.path.join(output_folder, str(nb_acq), 'AMICO',
                                  self.model)
            # Load dicom headers
            if not os.path.isfile(self.inputs['dcm']):
                err = "DICOM File %s not found."
                raise Exception(err % self.inputs['dcm'])
            sour_obj = dicom.read_file(self.inputs['dcm'])

            # Convert all niftis to dicoms
            convert_niftis_2_dicoms(
                outdir,
                sour_obj,
                osirix_folder,
                nb_acq)

            # Subtract the Cobj to the maps:
            subtract_obj_to_map(outdir, sour_obj, osirix_folder, nb_acq)

        # Make pdf:
        self.make_pdf()

        # Zip the DICOMs output:
        initdir = os.getcwd()
        # Zip all the files in the directory
        zip_name = os.path.join(self.jobdir, 'outputs', 'OsiriX', 'osirix.zip')
        os.chdir(os.path.join(self.jobdir, 'outputs', 'OsiriX'))
        os.system('zip -r %s * > /dev/null' % zip_name)
        # return to the initial directory:
        os.chdir(initdir)

        # Gzip nii:
        XnatUtils.gzip_nii(outdir)

    def make_pdf(self):
        """Method to make the PDF for the spider.

        :return: None
        """
        output_folder = os.path.join(self.jobdir, 'outputs')
        pdfs_dir = XnatUtils.makedir(os.path.join(output_folder, 'pdfs'))
        fpages = list()
        # Run matlab function
        for nb_acq in range(1, self.nb_acquisition + 1):
            pdf_page = os.path.join(output_folder, str(nb_acq),
                                    'VerdictMapAcq%d.pdf' % nb_acq)
            mat_lines = DEFAULT_PDF_MAKER.format(
                matlab_code=self.matlab_code,
                maps_folder=os.path.join(output_folder, str(nb_acq),
                                         'AMICO', self.model),
                subject=self.xnat_subject,
                output_folder=pdfs_dir,
                acq=nb_acq)
            matlab_script = os.path.join(output_folder,
                                         'run_pdf_page_%d.m' % nb_acq)
            with open(matlab_script, "w") as f:
                f.writelines(mat_lines)
            XnatUtils.run_matlab(matlab_script, verbose=True)
            # Get all PDFs:
            pdf_pages = XnatUtils.find_files(pdfs_dir, '.pdf')
            # Merge all pdfs into one:
            self.merge_pdf_pages(pdf_pages, pdf_page)
            fpages.append(pdf_page)

        if len(fpages) > 1:
            self.merge_pdf_pages(fpages, self.pdf_final)
        else:
            shutil.move(fpages[0], self.pdf_final)

    def finish(self):
        """Method to copy the results in dax.RESULTS_DIR."""
        results_dict = {'PDF': self.pdf_final,
                        'OsiriX': os.path.join(self.jobdir, 'outputs',
                                               'OsiriX', 'osirix.zip')}
        for nb_acq in range(1, self.nb_acquisition + 1):
            res = 'RMAPS%d' % nb_acq
            results_dict[res] = os.path.join(self.jobdir, 'outputs',
                                             str(nb_acq), 'AMICO', self.model)
        self.upload_dict(results_dict)
        self.end()

    def run_matlab(self, matlab_script, verbose=False):
        """Call MATLAB with -nodesktop -nosplash and -singlecompthread.

        :param matlab_script: Full path to the .m file to run
        :param verbose: True to print all MATLAB output to terminal, False to
         suppress.
        :return: None
        """
        self.time_writer("Matlab script: %s running ..." % matlab_script)
        cmd = "matlab -nodisplay -nodesktop -nojvm -nosplash -singleCompThread \
    < %s" % matlab_script
        if not verbose:
            matlabdir = os.path.dirname(matlab_script)
            prefix = os.path.basename(matlab_script).split('.')[0]
            cmd = '%s > %s' % (cmd, os.path.join(matlabdir, '%s_outlog.log' % prefix))
        os.system(cmd)
        self.time_writer("Matlab script: %s done" % matlab_script)


def write_dicom(pixel_array, filename, ds_ori,
                series_number, sop_id, series_description):
    """Write a dicom from a pixel_array (numpy).

    :param pixel_array: 2D numpy ndarray.
                        If pixel_array is larger than 2D, errors.
    :param filename: string name for the output file.
    :param ds_ori: original pydicom object of the pixel_array
    :param series_number: number of the series being processed
    :param sop_id: SOPInstanceUID for the DICOM
    :param series_description: series description for Osirix display
    """
    # Set the DICOM dataset
    file_meta = Dataset()
    file_meta.MediaStorageSOPClassUID = 'Secondary Capture Image Storage'
    file_meta.MediaStorageSOPInstanceUID = ds_ori.SOPInstanceUID
    file_meta.ImplementationClassUID = ds_ori.SOPClassUID
    ds = FileDataset(filename, {}, file_meta=file_meta, preamble="\0"*128)

    # Copy the tag from the original DICOM
    for tag, d_obj in ds_ori.items():
        if tag != ds_ori.data_element("PixelData").tag:
            ds[tag] = d_obj

    # Other tags to set
    ds.SeriesNumber = series_number
    sop_uid = sop_id + str(datetime.datetime.now()).replace('-', '')\
                                                   .replace(':', '')\
                                                   .replace('.', '')\
                                                   .replace(' ', '')
    ds.SOPInstanceUID = sop_uid[:-1]
    ds.ProtocolName = '%s Verdict MAP' % series_description
    # Set SeriesDate/ContentDate
    now = datetime.date.today()
    ds.SeriesDate = '%d%02d%02d' % (now.year, now.month, now.day)
    ds.ContentDate = '%d%02d%02d' % (now.year, now.month, now.day)
    ds.Modality = 'MR'
    ds.ConversionType = 'WSD'
    ds.StudyDescription = 'INNOVATE'
    ds.SeriesDescription = series_description
    ds.AcquisitionNumber = 1
    ds.SamplesperPixel = 1
    ds.PhotometricInterpretation = 'MONOCHROME2'
    ds.SecondaryCaptureDeviceManufctur = 'Python 2.7.3'
    nb_frames = pixel_array.shape[2]*pixel_array.shape[3]
    ds.NumberOfFrames = nb_frames
    ds.PixelRepresentation = 0
    ds.HighBit = 15
    ds.BitsStored = 8
    ds.BitsAllocated = 8
    ds.SmallestImagePixelValue = pixel_array.min()
    ds.LargestImagePixelValue = pixel_array.max()
    ds.Columns = pixel_array.shape[0]
    ds.Rows = pixel_array.shape[1]

    # Fixing the sequence if the number of frames was less than the original
    # it happens if we remove the last volume for phillips data (mean)
    if ds_ori.NumberOfFrames > nb_frames:
        new_seq = Sequence()
        for i in xrange(0, ds_ori.NumberOfFrames):
            if i % 5 == 0:  # take one slice for each (14)
                new_seq.append(ds_ori[0x5200, 0x9230][i])
        ds[0x5200, 0x9230].value = new_seq

    # Organise the array:
    pixel_array2 = np.zeros((pixel_array.shape[0]*pixel_array.shape[2],
                             pixel_array.shape[1]))
    for i in range(pixel_array.shape[2]):
        pixel_array2[pixel_array.shape[0]*i:pixel_array.shape[0]*(i+1),
                     :] = pixel_array[:, :, i, 0]

    # Set the Image pixel array
    if pixel_array2.dtype != np.uint8:
        pixel_array2 = pixel_array2.astype(np.uint8)

    ds.PixelData = pixel_array2.tostring()
    # Save the image
    ds.save_as(filename)


def convert_niftis_2_dicoms(nifti_folder, sour_obj, output_folder, nbacq):
    """Convert 3D niftis into DICOM files.

    :param nifti_path: path to the nifti file
    :param sour_obj: pydicom object for the source dicom
    :param output_folder: folder where the DICOM files will be saved
    :param nbacq: Acquisition number
    :return: None
    """
    if not os.path.isdir(nifti_folder):
        raise Exception("NIFTY Folder %s not found." % nifti_folder)
    for maps_name in ORDER_MAPS:
        # File
        nii_map = os.path.join(nifti_folder, '%s.nii' % maps_name)
        if not os.path.isfile(nii_map):
            raise Exception("NIFTY file %s not found." % nii_map)
        # Naming
        label = os.path.basename(nii_map)[:-4]
        series_description = '%s_%d_original' % (label, nbacq)

        # Edit Niftys
        f_img = nib.load(nii_map)
        f_data = f_img.get_data()
        # Rotation 270
        f_data = np.rot90(f_data)
        f_data = np.rot90(f_data)
        f_data = np.rot90(f_data)
        # Normalizing data:
        dmin = C_RANGE[label]['min']
        dmax = C_RANGE[label]['max']
        f_data[f_data < dmin] = dmin
        f_data[f_data > dmax] = dmax
        # Scale numbers to uint8
        f_data = (f_data - dmin)*255.0/(dmax - dmin)

        # Make output_folder:
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        # Series Number and SOP UID
        series_number = 87000 + C_RANGE[label]['index'] + int(nbacq) - 1
        sop_id = sour_obj.SOPInstanceUID.split('.')
        sop_id = '.'.join(sop_id[:-1])+'.'

        # Write the dicom
        filename = os.path.join(output_folder, '%s.dcm' % series_description)
        write_dicom(f_data, filename, sour_obj,
                    series_number, sop_id, series_description)


def subtract_obj_to_map(nii_folder, sour_obj, output_folder, nbacq):
    if not os.path.isdir(nii_folder):
        raise Exception("NIFTY Folder %s not found." % nii_folder)

    nii_mapobj = os.path.join(nii_folder, 'FIT_FobjCamino.nii')
    if not os.path.isfile(nii_mapobj):
        raise Exception("NIFTY file OBJ %s not found." % nii_mapobj)

    # Load image from NIFTI
    f_img_mobj = nib.load(nii_mapobj)
    f_data_obj = f_img_mobj.get_data()
    # Rotation 270
    f_data_obj = np.rot90(f_data_obj)
    f_data_obj = np.rot90(f_data_obj)
    f_data_obj = np.rot90(f_data_obj)
    mask = f_data_obj > C_RANGE['FIT_FobjCamino']['max']

    for index, map_name in enumerate(['FIT_fIC', 'FIT_cellularity',
                                      'FIT_fEES', 'FIT_fVASC', 'FIT_R']):
        # File
        nii_map = os.path.join(nii_folder, '%s.nii' % map_name)
        if not os.path.isfile(nii_map):
            raise Exception("NIFTY file %s not found." % nii_map)
        # Load image from NIFTI
        f_img = nib.load(nii_map)
        f_data = f_img.get_data()
        # Rotation 270
        f_data = np.rot90(f_data)
        f_data = np.rot90(f_data)
        f_data = np.rot90(f_data)
        # Normalizing data:
        dmin = C_RANGE[map_name]['min']
        dmax = C_RANGE[map_name]['max']
        f_data[f_data < dmin] = dmin
        f_data[f_data > dmax] = dmax

        # Scale numbers to uint8
        f_data = (f_data - dmin)*255.0/(dmax - dmin)

        # Subtract by setting all value to nan:
        f_data[mask] = np.nan

        # Make output_folder:
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        # map_name:
        series_description = '%s_%d_subtracted' % (map_name, nbacq)

        # Series Number and SOP UID
        series_number = 87000 + index + 5
        sop_id = sour_obj.SOPInstanceUID.split('.')
        sop_id = '.'.join(sop_id[:-1])+'.'

        filename = os.path.join(output_folder, '%s.dcm' % series_description)
        write_dicom(f_data, filename, sour_obj,
                    series_number, sop_id, series_description)


if __name__ == '__main__':
    args = parse_args()
    # generate spider object:
    spider_obj = Spider_Verdict(spider_path=sys.argv[0],
                                jobdir=args.temp_dir,
                                xnat_project=args.proj_label,
                                xnat_subject=args.subj_label,
                                xnat_session=args.sess_label,
                                proctype=args.proctype,
                                nb_acquisition=args.nb_acquisition,
                                matlab_code=args.matlab_code,
                                camino=args.camino,
                                spams=args.spams,
                                scheme_filename=args.scheme_filename,
                                model=args.model,
                                xnat_host=args.host,
                                xnat_user=args.user,
                                xnat_pass=None,
                                suffix=args.suffix)

    # print some information before starting
    spider_obj.print_init(args, "Benjamin Yvernault", "b.yvernault@ucl.ac.uk")

    # Pre-run method to download data from XNAT
    spider_obj.pre_run()

    # Run method
    spider_obj.run()

    # Finish method to copy results
    if not args.skipfinish:
        spider_obj.finish()
