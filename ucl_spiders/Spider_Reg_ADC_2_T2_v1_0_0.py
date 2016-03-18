'''
    Author:         Benjamin Yvernault
    contact:        b.yvernault@ucl.ac.uk
    Spider name:    Reg_ADC_2_T2
    Spider version: 1.0.0
    Creation date:  2016-03-15 13:33:03.911277
    Purpose:        Register ADC scan to T2 scan
'''

__author__ = "Benjamin Yvernault"
__email__ = "b.yvernault@ucl.ac.uk"
__purpose__ = "Register ADC scan to T2 scan"
__spider_name__ = "Reg_ADC_2_T2"
__version__ = "1.0.0"
__modifications__ = "2016-03-15 13:33:03.911277 - Original write"

# Python packages import
import os
import sys
import time
import dicom
import datetime
import numpy as np
import nibabel as nib
import subprocess as sb
import matplotlib.pyplot as plt
from dicom.dataset import Dataset, FileDataset
from dax import spiders, XnatUtils, SessionSpider

REG_ALADIN_CMD = "{exe_path} -ref {ref} -flo {flo} -res {res} -aff {aff} {args}"
REG_F3D_CMD = "{exe_path} -ref {ref} -flo {flo} -aff {aff} -cpp {cpp} -res {res} {args}"
DEFAULT_ARGS_REG_ALADIN = " -maxit 15 -ln 4 -lp 4 -interp 1"
DEFAULT_ARGS_REG_F3D = " -ln 4 -lp 4 -jl 0.1 -be 0.05 -maxit 250 -lncc 0 5.0 -sx 2.5"

def parse_args():
    '''
        Argument parser for the spider input variables
        by default (set in get_session_argparser):
            -p       : proj_label
            -s       : subj_label
            -e       : sess_label
            -d       : temp_dir
            --suffix : suffix (suffix for assessor proctype on XNAT)
            --host : host for XNAT (default: XNAT_HOST env variable)
            --user : user on XNAT (default: XNAT_USER env variable)
        your arguments:
            --adc : ADC scan ID in the session on XNAT
            --t2 : T2 scan ID in the session on XNAT
            --regAladin : path to reg_aladin's executable
            --argsRegAladin : arguments for Reg Aladin. Default: -maxit 15 -ln 4 -lp 4 -interp 1
            --regf3d : path to reg_f3d's executable
            --argRegf3d : arguments for reg_f3d. Default: -ln 4 -lp 4 -jl 0.1 -be 0.05 -maxit 250 -lncc 0 5.0 -sx 2.5
            --openmp_core : number of core use by reg_aladin. Default: one

        :return: argument parser object created by parse_args()
    '''
    ap = spiders.get_session_argparser("Reg_ADC_2_T2", "Register ADC scan to T2 scan")
    ap.add_argument("--adc", dest="adc_id", help="ADC Scan ID from XNAT.", required=True)
    ap.add_argument("--t2", dest="t2_id", help="T2 Scan ID from XNAT.", required=True)
    ap.add_argument("--regAladin", dest="reg_aladin_exe", help="path to reg_aladin's executable.", required=True)
    ap.add_argument("--argRegAladin", dest="args_reg_aladin", help="Argument for reg_aladin. Default: -maxit 15 -ln 4 -lp 4 -interp 1.", default=DEFAULT_ARGS_REG_ALADIN)
    ap.add_argument("--regf3d", dest="regf3d_exe", help="path to reg_f3d's executable.", required=True)
    ap.add_argument("--argRegf3d", dest="args_regf3d", help="Argument for reg_f3d. Default: -maxit 15 -ln 4 -lp 4 -interp 1.", default=DEFAULT_ARGS_REG_F3D)
    ap.add_argument("--openmp_core", dest="openmp_core", help="Number of core used by reg_aladin.", required=False, default=1)
    return ap.parse_args()

class Spider_Reg_ADC_2_T2(SessionSpider):
    '''
        Session Spider class to do: Register ADC scan to T2 scan

        :param spider_path: spider file path
        :param jobdir: directory for temporary files
        :param xnat_project: project ID on XNAT
        :param xnat_subject: subject label on XNAT
        :param xnat_session: experiment label on XNAT
        :param xnat_host: host for XNAT if not set in environment variables
        :param xnat_user: user for XNAT if not set in environment variables
        :param xnat_pass: password for XNAT if not set in environment variables
        :param suffix: suffix to the assessor creation
    '''
    def __init__(self, spider_path, jobdir, xnat_project, xnat_subject, xnat_session,
                 xnat_host=None, xnat_user=None, xnat_pass=None, suffix=""):
        super(Spider_Reg_ADC_2_T2, self).__init__(spider_path, jobdir, xnat_project, xnat_subject, xnat_session,
                                            xnat_host, xnat_user, xnat_pass, suffix)
        # Inputs
        self.adc_nii = list()
        self.t2_nii = list()
        self.adc_dcm = list()
        self.t2_dcm = list()

        # Outputs
        self.outputs = list()
        self.pdf_final = os.path.join(self.jobdir, 'Registration_ADC_2_T2.pdf')
        # Check Executable:
        self.reg_aladin_exe = self.check_exe(ARGS.reg_aladin_exe, 'reg_aladin')
        self.reg_f3d_exe = self.check_exe(ARGS.regf3d_exe, 'reg_f3d')
        # Print version for Niftyreg - GIFi
        pversion = sb.Popen([self.reg_aladin_exe, '--version'], stdout=sb.PIPE,
                                                                stderr=sb.PIPE)
        nve_version, _ = pversion.communicate()
        self.time_writer('reg_aladin (Niftyreg) version: %s' % (nve_version.strip()))
        # Print version for Niftyreg - GIFi
        pversion = sb.Popen([self.reg_f3d_exe, '--version'], stdout=sb.PIPE,
                                                                stderr=sb.PIPE)
        nve_version, _ = pversion.communicate()
        self.time_writer('reg_f3d (Niftyreg) version: %s' % (nve_version.strip()))

    def pre_run(self):
        '''
            Method to download data from XNAT

            :param argument_parse: argument parser object return by parse_args()
        '''
        resource = 'NIFTI' #resource to download from the scan on XNAT
        input_folder = os.path.join(self.jobdir, 'inputs')
        if not os.path.exists(input_folder):
            os.makedirs(input_folder)
        adc_folder = os.path.join(input_folder, 'ADC')
        if not os.path.exists(adc_folder):
            os.makedirs(adc_folder)
        t2_folder = os.path.join(input_folder, 'T2')
        if not os.path.exists(t2_folder):
            os.makedirs(t2_folder)
        adc_dcm = os.path.join(adc_folder, 'DICOM')
        if not os.path.exists(adc_dcm):
            os.makedirs(adc_dcm)
        t2_dcm = os.path.join(t2_folder, 'DICOM')
        if not os.path.exists(t2_dcm):
            os.makedirs(t2_dcm)

        self.adc_nii.extend(self.download(ARGS.adc_id, resource, adc_folder))
        self.t2_nii.extend(self.download(ARGS.t2_id, resource, t2_folder))
        self.adc_dcm.extend(self.download(ARGS.adc_id, 'DICOM', adc_dcm))
        self.t2_dcm.extend(self.download(ARGS.t2_id, 'DICOM', t2_dcm))


    @staticmethod
    def check_exe(executable, name):
        """
            Method to check the executable

            :param executable: executable path
            :param name: name of Executable
            :return: Complete path to the executable
        """
        if executable == name:
            print name
            return name
        else:
            executable = os.path.abspath(executable)
            if executable.endswith(name):
                path = executable
            elif os.path.isdir(executable):
                path = os.path.join(executable, name)

            if not os.path.exists(executable):
                raise Exception("Executable '%s' not found" % (executable))

            return path

    def run(self):
        '''
            Method running the process for the spider on the inputs data
        '''
        output_folder = os.path.join(self.jobdir, 'outputs')
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        # REG ALADIN:
        self.time_writer("reg_aladin with ref %s and flo %s" % (self.adc_nii[0], self.t2_nii[0]))
        aladin_output = os.path.join(output_folder, "ADC_reg_2_t2_reg_aladin.nii")
        affine_fpath = os.path.join(output_folder, "ADC_reg_2_t2_affine_transformation.txt")

        cmd = REG_ALADIN_CMD.format(exe_path=self.reg_aladin_exe,
                                    ref=self.t2_nii[0],
                                    flo=self.adc_nii[0],
                                    res=aladin_output,
                                    aff=affine_fpath,
                                    args=ARGS.args_reg_aladin)
        self.run_system_cmd(cmd)
        # Check that the affine file exists:
        if not os.path.exists(affine_fpath):
            raise Exception('Reg_aladin failed. File %s not found.' % affine_fpath)

        #REG_F3D
        self.time_writer("reg_f3d with ref %s and flo %s and aff %s" % (self.adc_nii[0], self.t2_nii[0], affine_fpath))
        f3d_output = os.path.join(output_folder, "ADC_reg_2_t2_reg_f3d.nii")
        f3d_cpp = os.path.join(output_folder, "ADC_reg_2_t2_reg_f3d_cpp.nii")
        cmd = REG_F3D_CMD.format(exe_path=self.reg_f3d_exe,
                                 ref=self.t2_nii[0],
                                 flo=self.adc_nii[0],
                                 res=f3d_output,
                                 cpp=f3d_cpp,
                                 aff=affine_fpath,
                                 args=ARGS.args_regf3d)
        self.run_system_cmd(cmd)
        XnatUtils.gzip_nii(output_folder)
        self.outputs.extend([{'label':'reg_aladin_results', 'image':aladin_output+'.gz'},
                             {'label':'reg_f3d_results', 'image':f3d_output+'.gz'}])
        # Make PDF
        self.make_pdf()
        # Generate DICOM version of the reg_f3d results:
        convert_nifti_2_dicoms(os.path.join(self.jobdir, 'outputs', "ADC_reg_2_t2_reg_f3d.nii.gz"),
                               self.t2_dcm[0], self.adc_dcm[0],
                               os.path.join(output_folder, 'reg_f3d_dicoms'), label="ADC_reg_f3d")

    def finish(self):
        '''
            Method to copy the results in the Spider Results folder dax.RESULTS_DIR
        '''
        results_dict = {'PDF': self.pdf_final,
                        'REG_ALA': os.path.join(self.jobdir, 'outputs', "ADC_reg_2_t2_reg_aladin.nii.gz"),
                        'AFF': os.path.join(self.jobdir, 'outputs', "ADC_reg_2_t2_affine_transformation.txt"),
                        'REG_F3D': os.path.join(self.jobdir, 'outputs', "ADC_reg_2_t2_reg_f3d.nii.gz"),
                        'OSIRIX' : os.path.join(self.jobdir, 'outputs', 'reg_f3d_dicoms'),
                        'CPP': os.path.join(self.jobdir, 'outputs', "ADC_reg_2_t2_reg_f3d_cpp.nii.gz")}
        self.upload_dict(results_dict)
        self.end()

    def plot_images(self, fig, data, number_images, subplot_index, label):
        '''
            Method to plot the images on the PDF for the first page

            :param fig: figure from matplotlib
            :param data: numpy array of the images (3D) to display
            :param number_images: number of images for the subplot (6)
            :param subplot_index: index of the line to display (1-6)
            :param label: y-label
            :return: None
        '''
        ax = fig.add_subplot(number_images, 3, 3*subplot_index+1)
        ax.imshow(np.rot90(data[:, :, data.shape[2]/2]), cmap = 'gray')
        if subplot_index == 0:
            ax.set_title('Axial', fontsize=7)
        ax.set_ylabel(label, fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
        ax = fig.add_subplot(number_images, 3, 3*subplot_index+2)
        ax.imshow(np.rot90(data[:, data.shape[1]/2, :]), cmap = 'gray')
        if subplot_index == 0:
            ax.set_title('Coronal', fontsize=7)
        ax.set_axis_off()
        ax = fig.add_subplot(number_images, 3, 3*subplot_index+3)
        ax.imshow(np.rot90(data[data.shape[0]/2, :, :]), cmap = 'gray')
        if subplot_index == 0:
            ax.set_title('Sagittal', fontsize=7)
        ax.set_axis_off()

    @staticmethod
    def save_pdf_page(fig, pdf_path, date, pdf_page_number, pdf_pages):
        """
            Method to save the pdf page to a file

            :param fig: figure from matplotlib
            :param pdf_path: numpy array of the images (3D) to display
            :param number_images: number of images for the subplot (6)
            :param subplot_index: index of the line to display (1-6)
            :param label: y-label
            :return: None
        """
        # Set footer and title
        fig.tight_layout()
        plt.figtext(0.5, 0.985, '-- Registration2Ref Pipeline PDF report --',
                    horizontalalignment='center', fontsize=10)
        footer = 'Date: %s -- page %s/%s -- PDF generated by TIG laboratory at UCL, London' % (str(date), str(pdf_page_number), str(pdf_pages))
        plt.figtext(0.5, 0.02, footer, horizontalalignment='center', fontsize=8)
        # Writing Figure
        #plt.show()
        fig.savefig(pdf_path, transparent=True, orientation='portrait', dpi=100)
        plt.close(fig)

    def make_pdf(self):
        '''
            Method to make the PDF for the spider

            :return: None
        '''
        ## Define output files
        # Variables:
        date = datetime.datetime.now()
        nb_images_per_page = 5

        # PDF path:
        pdf_pages_list = list()

        pdf_page_number = 0 # 5 images per pages
        fig = plt.figure(0, figsize=(7.5, 10))
        total_nb_pages = (len(self.outputs) + 2)/nb_images_per_page
        if total_nb_pages == 0:
            total_nb_pages = 1
        # Plot the Inputs
        # Open niftis with nibabel
        f_img = nib.load(self.t2_nii[0])
        f_img_data = f_img.get_data()
        # Draw
        if len(f_img_data.shape) == 3:
            data = f_img_data
        elif len(f_img_data.shape) == 4:
            data = f_img_data[:, :, :, f_img_data.shape[3]/2]
        self.plot_images(fig, data, nb_images_per_page, 0, 'T2 scan')
        # Open niftis with nibabel
        f_img = nib.load(self.adc_nii[0])
        f_img_data = f_img.get_data()
        # Draw
        if len(f_img_data.shape) == 3:
            data = f_img_data
        elif len(f_img_data.shape) == 4:
            data = f_img_data[:, :, :, f_img_data.shape[3]/2]
        self.plot_images(fig, data, nb_images_per_page, 1, 'ADC scan')

        page_0_saved = False
        for index, out_dict in enumerate(self.outputs):
            if not os.path.exists(out_dict['image']):
                raise Exception('%s output image not found.' % (out_dict['image']))

            if (index+2)/5 != pdf_page_number:
                page_0_saved = True
                pdf_name = 'Reg_ADC_2_T2_page%s.pdf' % (pdf_page_number)
                pdf_page = os.path.join(self.jobdir, pdf_name)
                # Save figure to PDF page
                self.save_pdf_page(fig, pdf_page, date, pdf_page_number, total_nb_pages)
                pdf_pages_list.append(pdf_page)
                pdf_page_number += 1
                fig = plt.figure(pdf_page_number, figsize=(7.5, 10))

            # Open niftis with nibabel
            f_img = nib.load(out_dict['image'])
            f_img_data = f_img.get_data()
            # Draw
            if len(f_img_data.shape) == 3:
                data = f_img_data
            elif len(f_img_data.shape) == 4:
                data = f_img_data[:, :, :, f_img_data.shape[3]/2]
            self.plot_images(fig, data, nb_images_per_page, index+2, out_dict['label'])

        # If only 5 or less images, save PDF zero
        if not page_0_saved:
            self.save_pdf_page(fig, self.pdf_final, date, 1, total_nb_pages)
        else:
            ## Join the two pages for the PDF:
            cmd = 'gs -q -sPAPERSIZE=letter -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=%s' % (self.pdf_final)
            for page in pdf_pages_list:
                cmd = '%s %s' % (cmd, page)
            self.time_writer('INFO:saving final PDF: %s ' % cmd)
            os.system(cmd)

def write_dicom(pixel_array, filename, ds_copy, ds_ori, volume_number,
                series_number, sop_id):
    """
    INPUTS:
    pixel_array: 2D numpy ndarray.  If pixel_array is larger than 2D, errors.
    filename: string name for the output file.
    """
    # Set to zero negatives values in the image:
    pixel_array[pixel_array<0] = 0

    # Set the DICOM dataset
    file_meta = Dataset()
    file_meta.MediaStorageSOPClassUID = 'Secondary Capture Image Storage'
    file_meta.MediaStorageSOPInstanceUID = ds_ori.SOPInstanceUID
    file_meta.ImplementationClassUID = ds_ori.SOPClassUID
    ds = FileDataset(filename, {}, file_meta = file_meta, preamble="\0"*128)

    # Copy the tag from the original DICOM
    for tag, value in ds_ori.items():
        if tag != ds_copy.data_element("PixelData").tag:
            ds[tag] = value

    # Other tags to set
    ds.SeriesNumber = series_number
    ds.SeriesDescription = ds_ori.SeriesDescription + ' reg_f3d'
    sop_uid = sop_id + str(datetime.datetime.now()).replace('-','').replace(':','').replace('.','').replace(' ','')
    ds.SOPInstanceUID = sop_uid[:-1]
    ds.ProtocolName = ds_ori.ProtocolName
    ds.InstanceNumber = volume_number
    # Set time:
    ti = time.time()
    ds.ContentDate = str(datetime.date.today()).replace('-','')
    ds.ContentTime = str(ti)

    # Copy from T2 the orientation tags:
    ds.PatientPosition = ds_copy.PatientPosition
    ds[0x20,0x32] = ds_copy[0x20,0x32] # Image Position
    ds[0x20,0x37] = ds_copy[0x20,0x37] # Image Orientation
    ds[0x18,0x50] = ds_copy[0x18,0x50] # Slice Thicknes
    ds[0x18,0x88] = ds_copy[0x18,0x88] # Spacing Between Slices
    ds[0x18,0x1312] = ds_copy[0x18,0x1312] # In-plane Phase Encoding
    ds[0x28,0x10] = ds_copy[0x28,0x10] # rows
    ds[0x28,0x11] = ds_copy[0x28,0x11] # columns
    ds[0x28,0x30] = ds_copy[0x28,0x30] # Pixel spacing

    # Set the Image pixel array
    if pixel_array.dtype != np.uint16:
        pixel_array = pixel_array.astype(np.uint16)
    ds.PixelData = pixel_array.tostring()

    # Save the image
    ds.save_as(filename)

def convert_nifti_2_dicoms(nifti_path, dicom_target, dicom_source, output_folder, label=None):
    """
    Convert 4D niftis generated by reg_f3d into DICOM files.

    :param nifti_path: path to the nifti file
    :param dicom_source: one dicom file from the source
     for the registration for header info
    :param dicom_target: one dicom file from the target
     for the registration for header info
    :param output_folder: folder where the DICOM files will be saved
    :param label: name for the output dicom files
    :return: None
    """
    if not os.path.isfile(nifti_path):
        raise Exception("File %s not found after reg_f3d." % nifti_path)
    # Load image from NIFTI
    f_img = nib.load(nifti_path)
    f_img_data = f_img.get_data()

    # Load dicom headers
    if not os.path.isfile(dicom_target):
        raise Exception("DICOM File %s not found after reg_f3d." % dicom_target)
    t2_dcm_obj = dicom.read_file(dicom_target)
    if not os.path.isfile(dicom_source):
        raise Exception("DICOM File %s not found after reg_f3d." % dicom_source)
    adc_dcm_obj = dicom.read_file(dicom_source)

    # Make output_folder:
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Series Number and SOP UID
    ti = time.time()
    series_number = 86532 + int(str(ti)[2:4]) + int(str(ti)[4:6])
    sop_id = adc_dcm_obj.SOPInstanceUID.split('.')
    sop_id = '.'.join(sop_id[:-1])+'.'

    for volume_index in range(f_img_data.shape[2]):
        filename = os.path.join(output_folder, '%s_%s.dcm' % (label, str(volume_index)))
        write_dicom(np.rot90(f_img_data[:, :, volume_index]), filename,
                    t2_dcm_obj, adc_dcm_obj, volume_index, series_number, sop_id)

if __name__ == '__main__':
    ARGS = parse_args()
    # generate spider object:
    spider_obj = Spider_Reg_ADC_2_T2(spider_path=sys.argv[0],
                                     jobdir=ARGS.temp_dir,
                                     xnat_project=ARGS.proj_label,
                                     xnat_subject=ARGS.subj_label,
                                     xnat_session=ARGS.sess_label,
                                     xnat_host=ARGS.host,
                                     xnat_user=ARGS.user,
                                     xnat_pass=None,
                                     suffix=ARGS.suffix)

    # print some information before starting
    spider_obj.print_init(ARGS, "Benjamin Yvernault", "b.yvernault@ucl.ac.uk")

    # Pre-run method to download data from XNAT
    spider_obj.pre_run()

    # Run method
    spider_obj.run()

    # Finish method to copy results
    spider_obj.finish()
