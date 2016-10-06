"""Spider_Comput_ADC_Verdict.

Author:         Benjamin Yvernault
contact:        b.yvernault@ucl.ac.uk
Spider name:    Verdict
Spider version: 1.0.0
Creation date:  2016-08-10 17:48:51.701453
Purpose:        Generate ADC from Verdict files after registration.
"""

# Python packages import
import os
import sys
import dicom
import datetime
import numpy as np
import nibabel as nib
from dicom.sequence import Sequence
from dicom.dataset import Dataset, FileDataset
from dax import XnatUtils, spiders, SessionSpider

__author__ = "Benjamin Yvernault"
__email__ = "b.yvernault@ucl.ac.uk"
__purpose__ = "Generate ADC from Verdict files after registration."
__spider_name__ = "Verdict"
__version__ = "1.0.0"
__modifications__ = """2016-08-10 17:48:51.701453 - Original write"""


DEFAULT_SCHEME_FILE = "NOptimisedADC_IN.scheme"
DEFAULT_VERDICT_TEMPLATE = """
addpath(genpath('{matlab_code}'));
compute_ADC_VERDICT('{input_path}',\
'{output}',\
'{scheme_filename}',\
'{camino}');
"""
DICOM_SCAN_TYPE = ['WIP b3000_90 SENSE', 'SWITCH DB TO YES b3000_80',
                   'b3000_80']


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
    ap = spiders.get_session_argparser("Compute_ADC_Verdict", __purpose__)
    ap.add_argument("--proctype", dest="proctype", default=None, required=True,
                    help="Assessor type containing the registered nifti.")
    ap.add_argument("--nbAcq", dest="nb_acquisition", default=1,
                    required=True,
                    help="Number of Acquisition of VERDICT scans (1 or 2).")
    ap.add_argument("--mc", dest="matlab_code", default=None, required=True,
                    help="Matlab code folder where is \
launch_AMICO_for_INNOVATE.")
    ap.add_argument("--camino", dest="camino", default=None, required=True,
                    help="Path to Camino folder.")
    ap.add_argument("--scheme", dest="scheme_filename", default=None,
                    help="Path to schemeFilename \
(NOptimisedV_IN.scheme).")
    return ap.parse_args()


class Spider_Comput_ADC_Verdict(SessionSpider):
    """Session Spider: Spider_Comput_ADC_Verdict.

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
                 matlab_code, camino, scheme_filename=None,
                 xnat_host=None, xnat_user=None, xnat_pass=None,
                 suffix=""):
        """Entry point for Spider_Verdict Class."""
        super(Spider_Comput_ADC_Verdict,
              self).__init__(spider_path, jobdir,
                             xnat_project, xnat_subject, xnat_session,
                             xnat_host, xnat_user, xnat_pass,
                             suffix)
        self.proctype = proctype
        self.nb_acquisition = int(nb_acquisition)

        self.inputs = dict()
        self.pdf_final = os.path.join(self.jobdir,
                                      '%s_Compute_ADC_VERDICT_report.pdf'
                                      % xnat_session)

        self.matlab_code = matlab_code
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
                    input_path=self.inputs[nb_acq],
                    output=folder,
                    scheme_filename=self.scheme_filename,
                    camino=self.camino)
            matlab_script = os.path.join(folder, 'run_matlab_verdict.m')
            with open(matlab_script, "w") as f:
                f.writelines(mat_lines)
            self.run_matlab(matlab_script, verbose=True)

            # Generate Dicom for OsiriX
            out_nii = os.path.join(folder, 'FIT_ADC.nii')

            # Load dicom headers
            if not os.path.isfile(self.inputs['dcm']):
                err = "DICOM File %s not found."
                raise Exception(err % self.inputs['dcm'])
            sour_obj = dicom.read_file(self.inputs['dcm'])

            # Convert all niftis to dicoms
            convert_nifti_2_dicoms(
                out_nii,
                sour_obj,
                osirix_folder,
                nb_acq)

            # Gzip nii:
            XnatUtils.gzip_nii(folder)

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

    def make_pdf(self):
        """Method to make the PDF for the spider.

        :return: None
        """
        list_slices = [3, 6, 9, 12]
        slices = {'0': list_slices}
        labels = {'0': 'ADC Acquisition 1'}
        vmins = {'0': 0}
        vmaxs = {'0': 5*10**-9}
        images = [os.path.join(self.jobdir, 'outputs', '1', 'FIT_ADC.nii.gz')]
        if self.nb_acquisition == 2:
            images.append(os.path.join(self.jobdir, 'outputs', '2',
                                       'FIT_ADC.nii.gz'))
            labels['1'] = 'ADC Acquisition 2'
            slices['1'] = list_slices
            vmins['1'] = 0
            vmaxs['1'] = 5*10**-9
        self.plot_images_page(self.pdf_final, 1, images, 'Compute ADC Verdict',
                              image_labels=labels, slices=slices,
                              vmins=vmins, vmaxs=vmaxs)

    def finish(self):
        """Method to copy the results in dax.RESULTS_DIR."""
        results_dict = {'PDF': self.pdf_final,
                        'OsiriX': os.path.join(self.jobdir, 'outputs',
                                               'OsiriX', 'osirix.zip')}
        for nb_acq in range(1, self.nb_acquisition+1):
            acq_folder = os.path.join(self.jobdir, 'outputs', str(nb_acq))
            res = 'ADC%d' % nb_acq
            results_dict[res] = os.path.join(acq_folder, 'FIT_ADC.nii.gz')
            files = [os.path.join(acq_folder, 'prepareADC.Bfloat'),
                     os.path.join(acq_folder, 'temporalADC.Bfloat'),
                     os.path.join(acq_folder, 'run_matlab_verdict.m')]
            results_dict['tmp_%d' % nb_acq] = files

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
            cmd = cmd+' > '+os.path.join(matlabdir, prefix+'_outlog.log')
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
    ds.ProtocolName = '%s Verdict' % series_description
    # Set SeriesDate/ContentDate
    now = datetime.date.today()
    ds.SeriesDate = '%d%02d%02d' % (now.year, now.month, now.day)
    ds.ContentDate = '%d%02d%02d' % (now.year, now.month, now.day)
    ds.Modality = 'MR'
    ds.StudyDescription = 'INNOVATE'
    ds.SeriesDescription = series_description
    ds.AcquisitionNumber = 1
    ds.SamplesperPixel = 1
    ds.PhotometricInterpretation = 'MONOCHROME2'
    ds.SecondaryCaptureDeviceManufctur = 'Python 2.7.11'
    ds.NumberOfFrames = pixel_array.shape[2]
    ds.PixelRepresentation = 0
    ds.HighBit = 15
    ds.BitsStored = 8
    ds.BitsAllocated = 8
    ds.SmallestImagePixelValue = pixel_array.min()
    ds.LargestImagePixelValue = pixel_array.max()
    ds.Columns = pixel_array.shape[0]
    ds.Rows = pixel_array.shape[1]

    # Organise the array:
    pixel_array2 = np.zeros((pixel_array.shape[0]*pixel_array.shape[2],
                             pixel_array.shape[1]))
    for i in range(pixel_array.shape[2]):
        pixel_array2[pixel_array.shape[0]*i:pixel_array.shape[0]*(i+1),
                     :] = pixel_array[:, :, i]

    # Set the Image pixel array
    if pixel_array2.dtype != np.uint8:
        pixel_array2 = pixel_array2.astype(np.uint8)

    ds.PixelData = pixel_array2.tostring()
    # Save the image
    ds.save_as(filename)


def convert_nifti_2_dicoms(nifti_file, sour_obj, output_folder, nbacq):
    """Convert 3D niftis into DICOM files.

    :param nifti_file: path to the nifti file
    :param sour_obj: pydicom object for the source dicom
    :param output_folder: folder where the DICOM files will be saved
    :param nbacq: Acquisition number
    :return: None
    """
    if not os.path.isfile(nifti_file):
        raise Exception("NIFTY file %s not found." % nifti_file)
    # Naming
    label = os.path.basename(nifti_file)[:-4]
    series_description = '%s_%d_original' % (label, nbacq)

    # Edit Niftys
    f_img = nib.load(nifti_file)
    f_data = f_img.get_data()
    # Rotation 270
    f_data = np.rot90(f_data)
    f_data = np.rot90(f_data)
    f_data = np.rot90(f_data)
    # Scale numbers to uint8
    # f_data = (f_data - f_data.min())*255.0/(f_data.max() - f_data.min())
    # Normalizing data:
    dmin = 0
    dmax = 5*10**-9
    f_data[f_data < dmin] = dmin
    f_data[f_data > dmax] = dmax

    # Scale numbers to uint8
    f_data = (f_data - dmin)*255.0/(dmax - dmin)

    # Subtract by setting all value to nan:
    # f_data[mask] = np.nan

    # Make output_folder:
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Series Number and SOP UID
    series_number = 88000 + nbacq
    sop_id = sour_obj.SOPInstanceUID.split('.')
    sop_id = '.'.join(sop_id[:-1])+'.'

    # Write the dicom
    filename = os.path.join(output_folder, '%s.dcm' % series_description)
    write_dicom(f_data, filename, sour_obj,
                series_number, sop_id, series_description)


if __name__ == '__main__':
    args = parse_args()
    # generate spider object:
    spider_obj = Spider_Comput_ADC_Verdict(
                                spider_path=sys.argv[0],
                                jobdir=args.temp_dir,
                                xnat_project=args.proj_label,
                                xnat_subject=args.subj_label,
                                xnat_session=args.sess_label,
                                proctype=args.proctype,
                                nb_acquisition=args.nb_acquisition,
                                matlab_code=args.matlab_code,
                                camino=args.camino,
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
