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
import time
import dicom
import shutil
import datetime
import numpy as np
import nibabel as nib
from dicom.dataset import Dataset, FileDataset
from dax import XnatUtils, spiders, SessionSpider

__author__ = "Benjamin Yvernault"
__email__ = "b.yvernault@ucl.ac.uk"
__purpose__ = "Generate Verdict Map from all Verdict scans registered together \
and merge as one nifti"
__spider_name__ = "Verdict"
__version__ = "1.0.0"
__modifications__ = """2016-08-10 17:48:51.701453 - Original write"""


DEFAULT_SCHEME_FILE = "NOptimisedV_IN.scheme"
DEFAULT_VERDICT_TEMPLATE = """
addpath(genpath('{matlab_code}'));
launch_AMICO_for_INNOVATE('{input_path}',\
'{subject}',\
'{filename}',\
'{output}',\
'{project}',\
'{amico}',\
'{camino}',\
'{spams}',\
'{scheme_filename}');
"""
DEFAULT_PDF_MAKER = """
% Code to generate PDF for VERDICT MAP

% Add path where matlab Code is
addpath(genpath('{matlab_code}'));

% Maps Name:
maps_name = {{'fIC','R','cellularity','fEES','fVASC','FobjCamino'}};

% Open an image to see the number of slices
f_file = fullfile('{maps_folder}', 'FIT_FobjCamino.nii');
aux = load_untouch_nii(f_file);
nb_slices = aux.hdr.dime.dim(4);
for i=1:nb_slices
    filename = ['DisplayMapsVerdictSlice' num2str(i, '%02d');];
    jpg_path = fullfile('{output_folder}', [filename '.pdf']);
    plot_oneslice_selectedmaps('{maps_folder}','{subject}',i,maps_name,1,1,1,jpg_path);
end
"""
DEFAULT_MAPS = ['FIT_FobjCamino.nii', 'FIT_R.nii', 'FIT_fEES.nii',
                'FIT_fIC.nii', 'FIT_fVASC.nii', 'FIT_cellularity.nii']
DICOM_SCAN_TYPE = ['SWITCH DB TO YES b3000_80']
MAX = {'FIT_FobjCamino': 50,
       'FIT_R': 1,
       'FIT_fEES': 50,
       'FIT_fIC': 1,
       'FIT_fVASC': 1,
       'FIT_cellularity': 50}


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
        --amico : path to AMICO folder
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
launch_AMICO_for_INNOVATE.")
    ap.add_argument("--amico", dest="amico", default=None, required=True,
                    help="Path to AMICO folder.")
    ap.add_argument("--camino", dest="camino", default=None, required=True,
                    help="Path to Camino folder.")
    ap.add_argument("--spams", dest="spams", default=None, required=True,
                    help="Path to spams-matlab folder.")
    ap.add_argument("--scheme", dest="scheme_filename", default=None,
                    help="Path to schemeFilename \
(NOptimisedV_IN.scheme).")
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
                 xnat_session, proctype, nb_acquisition,
                 matlab_code, amico, camino, spams, scheme_filename=None,
                 xnat_host=None, xnat_user=None, xnat_pass=None,
                 suffix=""):
        """Entry point for Spider_Verdict Class."""
        super(Spider_Verdict,
              self).__init__(spider_path, jobdir,
                             xnat_project, xnat_subject, xnat_session,
                             xnat_host, xnat_user, xnat_pass,
                             suffix)
        self.proctype = proctype
        self.nb_acquisition = int(nb_acquisition)

        self.inputs = dict()
        self.pdf_final = os.path.join(self.jobdir, 'VERDICT_report.pdf')

        self.matlab_code = matlab_code
        self.amico = amico
        self.spams = spams
        self.camino = camino
        if not scheme_filename:
            self.scheme_filename = os.path.join(matlab_code,
                                                DEFAULT_SCHEME_FILE)
        else:
            self.scheme_filename = scheme_filename

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
            file_gzip = XnatUtils.download_file_from_obj(
                                  folder, a.out_resource(res_name))
            XnatUtils.ungzip_nii(folder)
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
                    amico=self.amico,
                    camino=self.camino,
                    spams=self.spams,
                    scheme_filename=self.scheme_filename)
            matlab_script = os.path.join(output_folder, 'run_verdict.m')
            with open(matlab_script, "w") as f:
                f.writelines(mat_lines)
            XnatUtils.run_matlab(matlab_script)

            # Generate Dicom for OsiriX
            outdir = os.path.join(output_folder, str(nb_acq), 'AMICO',
                                  'VerdictProstate_Rmaps')
            # Load dicom headers
            if not os.path.isfile(self.inputs['dcm']):
                err = "DICOM File %s not found."
                raise Exception(err % self.inputs['dcm'])
            sour_obj = dicom.read_file(self.inputs['dcm'])
            for maps in DEFAULT_MAPS:
                nii = os.path.join(outdir, maps)
                convert_nifti_2_dicoms(
                    nii,
                    sour_obj,
                    osirix_folder,
                    os.path.basename(nii)[:-4],
                    label=os.path.basename(nii)[:-4])

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
        for nb_acq in range(1, self.nb_acquisition+1):
            pdf_page = 'VerdictMapAcq%d.pdf' % nb_acq
            mat_lines = DEFAULT_PDF_MAKER.format(
                    matlab_code=self.matlab_code,
                    maps_folder=os.path.join(output_folder, str(nb_acq),
                                             'AMICO',
                                             'VerdictProstate_Rmaps'),
                    subject=self.xnat_subject,
                    output_folder=pdfs_dir)
            matlab_script = os.path.join(output_folder,
                                         'run_pdf_page_%d.m' % nb_acq)
            with open(matlab_script, "w") as f:
                f.writelines(mat_lines)
            XnatUtils.run_matlab(matlab_script)
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
                                               'OsiriX')}
        for nb_acq in range(1, self.nb_acquisition+1):
            res = 'RMAPS%d' % nb_acq
            results_dict[res] = os.path.join(self.jobdir, 'outputs',
                                             str(nb_acq), 'AMICO',
                                             'VerdictProstate_Rmaps')
        self.upload_dict(results_dict)
        self.end()


def write_dicom(pixel_array, filename, ds_ori, volume_number,
                series_number, sop_id, series_description):
    """Write a dicom from a pixel_array (numpy).

    :param pixel_array: 2D numpy ndarray.
                        If pixel_array is larger than 2D, errors.
    :param filename: string name for the output file.
    :param ds_ori: original pydicom object of the pixel_array
    :param volume_number: number of the volume being processed
    :param series_number: number of the series being processed
    :param sop_id: SOPInstanceUID for the DICOM
    :param series_description: series description for Osirix display
    """
    # Set to zero negatives values in the image:
    # pixel_array[pixel_array < 0] = 0

    # Set the DICOM dataset
    file_meta = Dataset()
    file_meta.MediaStorageSOPClassUID = 'Secondary Capture Image Storage'
    file_meta.MediaStorageSOPInstanceUID = ds_ori.SOPInstanceUID
    file_meta.ImplementationClassUID = ds_ori.SOPClassUID
    ds = FileDataset(filename, {}, file_meta=file_meta, preamble="\0"*128)

    # Copy the tag from the original DICOM
    for tag, value in ds_ori.items():
        if tag != ds_ori.data_element("PixelData").tag:
            ds[tag] = value

    # Other tags to set
    ds.SeriesNumber = series_number
    ds.SeriesDescription = series_description
    sop_uid = sop_id + str(datetime.datetime.now()).replace('-', '')\
                                                   .replace(':', '')\
                                                   .replace('.', '')\
                                                   .replace(' ', '')
    ds.SOPInstanceUID = sop_uid[:-1]
    ds.ProtocolName = series_description
    ds.InstanceNumber = volume_number
    # Number of Frames:
    del ds[0x0028, 0x0008]

    # Copy rows/columns
    ds[0x00280010].value = pixel_array.shape[0]  # rows
    ds[0x00280011].value = pixel_array.shape[1]  # columns

    # Set the Image pixel array
    if pixel_array.dtype != np.uint16:
        pixel_array = pixel_array.astype(np.uint16)
    ds.PixelData = pixel_array.tostring()

    # Save the image
    ds.save_as(filename)


def convert_nifti_2_dicoms(nifti_path, sour_obj, output_folder,
                           series_description, label=None):
    """Convert 3D niftis into DICOM files.

    :param nifti_path: path to the nifti file
    :param sour_obj: pydicom object for the source dicom
    :param output_folder: folder where the DICOM files will be saved
    :param series_description: series description for the scan
    :param label: name for the output dicom files
    :return: None
    """
    if not os.path.isfile(nifti_path):
        raise Exception("NIFTY File %s not found." % nifti_path)
    # Load image from NIFTI
    f_img = nib.load(nifti_path)
    f_img_data = f_img.get_data()
    # Normalise the image for better display:
    minimun = 0
    f_img_data[f_img_data < minimun] = minimun
    f_img_data[f_img_data > MAX.get(label, 1.0)] = MAX.get(label, 1.0)
    f_img_data *= 255.0/np.float64(f_img_data.max())

    # Make output_folder:
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Series Number and SOP UID
    str_ti = "%f" % time.time()
    series_number = 86532 + int(str_ti[-4:-2]) + int(str_ti[-2:])
    sop_id = sour_obj.SOPInstanceUID.split('.')
    sop_id = '.'.join(sop_id[:-1])+'.'

    for slice_index in range(f_img_data.shape[2]):
        filename = os.path.join(output_folder, '%s_%02d.dcm' %
                                               (label, slice_index))
        write_dicom(np.rot90(f_img_data[:, :, slice_index]),
                    filename, sour_obj, slice_index,
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
                                amico=args.amico,
                                camino=args.camino,
                                spams=args.spams,
                                scheme_filename=args.scheme_filename,
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
    spider_obj.finish()
