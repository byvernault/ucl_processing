"""Spider to regirster different modalities for prostate cancer project.

Author:         Benjamin Yvernault
contact:        b.yvernault@ucl.ac.uk
Spider name:    Registration_Prostate
Spider version: 1.0.0
Creation date:  2016-03-15 13:33:03.911277
Purpose:        Register ADC scan to T2 scan
"""

# Python packages import
import os
import sys
import glob
import time
import dicom
import shutil
import datetime
import numpy as np
import nibabel as nib
import subprocess as sb
import matplotlib.pyplot as plt
from dicom.dataset import Dataset, FileDataset
from dax import spiders, XnatUtils, SessionSpider

__author__ = "Benjamin Yvernault"
__email__ = "b.yvernault@ucl.ac.uk"
__purpose__ = "Register Prostate Scans (DWI-ADC and DCE to T2)"
__spider_name__ = "Reg_ADC_2_T2"
__version__ = "1.0.0"
__modifications__ = """2016-03-15 13:33:03.911277 - Original write
2016-05-10 18:04:01 - Update to new format respecting pep8
"""

REG_ALADIN_CMD = "{exe_path} \
-ref {ref} \
-flo {flo} \
-res {res} \
-aff {aff} \
{args}"
REG_F3D_CMD = "{exe_path} \
-ref {ref} \
-flo {flo} \
-aff {aff} \
-cpp {cpp} \
-res {res} \
{args}"
DEFAULT_ARGS_REG_ALADIN = " -maxit 15 -ln 4 -lp 4 -interp 1"
DEFAULT_ARGS_REG_F3D = " -ln 4 -lp 4 -jl 0.1 \
-be 0.05 -maxit 250 -lncc 0 5.0 -sx 2.5"
# DICOMs TAG to copy
TAGS_TO_COPY = [0x00185100,  # Patient Position
                0x00180050,  # Slice Thicknes
                0x00180088,  # Spacing Between Slices
                0x00181312,  # In-plane Phase Encoding
                0x00200032,  # Image Position
                0x00200037,  # Image Orientation
                0x00201041,  # Slice Location
                0x00280010,  # rows
                0x00280011,  # columns
                0x00280030]  # Pixel spacing


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
        --adc : ADC scan ID in the session on XNAT
        --t2 : T2 scan ID in the session on XNAT
        --regAladin : path to reg_aladin's executable
        --argsRegAladin : arguments for Reg Aladin.
          Default: -maxit 15 -ln 4 -lp 4 -interp 1
        --regf3d : path to reg_f3d's executable
        --argRegf3d : arguments for reg_f3d.
          Default: -ln 4 -lp 4 -jl 0.1 -be 0.05 -maxit 250 -lncc 0 5.0 -sx 2.5
        --openmp_core : number of core use by reg_aladin. Default: one

    :return: argument parser object created by parse_args()
    """
    ap = spiders.get_session_argparser("Spider_Registration_Prostate",
                                       __purpose__)
    ap.add_argument("--sources", dest="sources_id", required=True,
                    help="Source Scans ID from XNAT.")
    ap.add_argument("--target", dest="target_id", required=True,
                    help="Target Scan ID from XNAT.")
    ap.add_argument("--regAladin", dest="reg_aladin_exe", required=True,
                    help="path to reg_aladin's executable.")
    ap.add_argument("--argsRegAladin", dest="args_reg_aladin",
                    help="Argument for reg_aladin. \
Default: -maxit 15 -ln 4 -lp 4 -interp 1.", default=DEFAULT_ARGS_REG_ALADIN)
    ap.add_argument("--regf3d", dest="regf3d_exe", required=True,
                    help="path to reg_f3d's executable.")
    ap.add_argument("--argsRegf3d", dest="args_regf3d",
                    help="Argument for reg_f3d. \
Default: -maxit 15 -ln 4 -lp 4 -interp 1.", default=DEFAULT_ARGS_REG_F3D)
    ap.add_argument("--openmp_core", dest="openmp_core", default=1,
                    help="Number of core used by reg_aladin.")
    return ap.parse_args()


class Spider_Registration_Prostate(SessionSpider):
    """Session Spider class to do: Register ADC scan to T2 scan.

    :param spider_path: spider file path
    :param jobdir: directory for temporary files
    :param xnat_project: project ID on XNAT
    :param xnat_subject: subject label on XNAT
    :param xnat_session: experiment label on XNAT
    :param xnat_host: host for XNAT if not set in environment variables
    :param xnat_user: user for XNAT if not set in environment variables
    :param xnat_pass: password for XNAT if not set in environment variables
    :param suffix: suffix to the assessor creation
    """

    def __init__(self, spider_path, jobdir,
                 xnat_project, xnat_subject, xnat_session,
                 xnat_host=None, xnat_user=None, xnat_pass=None, suffix=""):
        """Entry point for Spider_Registration_Prostate Class."""
        super(Spider_Registration_Prostate, self).__init__(spider_path, jobdir,
                                                           xnat_project,
                                                           xnat_subject,
                                                           xnat_session,
                                                           xnat_host,
                                                           xnat_user,
                                                           xnat_pass, suffix)
        # Inputs
        self.target = dict()
        self.sources = dict()

        # Outputs
        self.outputs = dict()
        self.pdf_final = os.path.join(self.jobdir, 'Registration_prostate.pdf')
        # Check Executable:
        self.reg_aladin_exe = self.check_exe(ARGS.reg_aladin_exe, 'reg_aladin')
        self.reg_f3d_exe = self.check_exe(ARGS.regf3d_exe, 'reg_f3d')
        # Print version for Niftyreg - GIFi
        pversion = sb.Popen([self.reg_aladin_exe, '--version'],
                            stdout=sb.PIPE,
                            stderr=sb.PIPE)
        nve_version, _ = pversion.communicate()
        self.time_writer('reg_aladin (Niftyreg) version: %s' %
                         (nve_version.strip()))
        # Print version for Niftyreg - GIFi
        pversion = sb.Popen([self.reg_f3d_exe, '--version'],
                            stdout=sb.PIPE,
                            stderr=sb.PIPE)
        nve_version, _ = pversion.communicate()
        self.time_writer('reg_f3d (Niftyreg) version: %s' %
                         (nve_version.strip()))

    def pre_run(self):
        """Method to download data from XNAT.

        :param argument_parse: argument parser object return by parse_args()
        """
        # Make directory
        input_folder = XnatUtils.makedir(os.path.join(self.jobdir, 'inputs'),
                                         subdir=False)

        # Target
        target_folder = XnatUtils.makedir(os.path.join(input_folder,
                                                       ARGS.target_id),
                                          subdir=False)
        target_dcm = XnatUtils.makedir(os.path.join(target_folder, 'DICOM'),
                                       subdir=False)
        self.time_writer('Connection to XNAT')
        xnat = XnatUtils.get_interface(host=self.host,
                                       user=self.user,
                                       pwd=self.pwd)
        self.time_writer('Downloading target %s ...' % ARGS.target_id)
        target_scan = XnatUtils.select_obj(xnat,
                                           ARGS.proj_label,
                                           ARGS.subj_label,
                                           ARGS.sess_label,
                                           ARGS.target_id)
        tnii_obj = target_scan.resource('NIFTI')
        self.target['nii'] = XnatUtils.download_file_from_obj(target_folder,
                                                              tnii_obj)
        tdcm_obj = target_scan.resource('DICOM')
        self.target['dcm'] = XnatUtils.download_files_from_obj(target_dcm,
                                                               tdcm_obj)
        self.target['type'] = target_scan.attrs.get('type')
        self.target['ID'] = ARGS.target_id

        # Sources
        sources_list = XnatUtils.get_input_list(ARGS.sources_id, list())
        self.time_writer('Downloading sources %s ...' % sources_list)
        for scan_id in sources_list:
            # Make directories
            spath = os.path.join(input_folder, scan_id)
            source_folder = XnatUtils.makedir(spath, subdir=False)
            dpath = os.path.join(source_folder, 'DICOM')
            source_dcm = XnatUtils.makedir(dpath, subdir=False)
            source_scan = XnatUtils.select_obj(xnat,
                                               ARGS.proj_label,
                                               ARGS.subj_label,
                                               ARGS.sess_label,
                                               scan_id)
            snii_obj = source_scan.resource('NIFTI')
            nii_list = XnatUtils.download_file_from_obj(source_folder,
                                                        snii_obj)
            sdcm_obj = source_scan.resource('DICOM')
            dcm_list = XnatUtils.download_file_from_obj(source_dcm, sdcm_obj)
            self.sources[scan_id] = dict()
            self.sources[scan_id]['nii'] = nii_list
            self.sources[scan_id]['dcm'] = dcm_list
            self.sources[scan_id]['type'] = source_scan.attrs.get('type')
            self.sources[scan_id]['ID'] = scan_id

        xnat.disconnect()
        self.time_writer('Disconnection of XNAT')

    @staticmethod
    def check_exe(executable, name):
        """Method to check the executable.

        :param executable: executable path
        :param name: name of Executable
        :return: Complete path to the executable
        """
        if executable == name:
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
        """Method running the process for the spider on the inputs data."""
        output_folder = XnatUtils.makedir(os.path.join(self.jobdir, 'outputs'),
                                          subdir=False)

        # Sort the DICOM T2 to convert the registered modalities NIFTI to DICOM
        dcm_obj_sorted = dict()
        for dcm_file in self.target['dcm']:
            # Load dicom headers
            if not os.path.isfile(dcm_file):
                err = "DICOM File %s not found after download."
                raise Exception(err % dcm_file)
            t2_dcm_obj = dicom.read_file(dcm_file)
            dcm_obj_sorted[t2_dcm_obj[0x00200032].value[2]] = t2_dcm_obj
        dcm_obj_sorted_list = [dcm_obj_sorted[key] for key in
                               sorted(dcm_obj_sorted)]

        # REG ALADIN:
        for scan_id, res_dict in self.sources.items():
            # Organise folders for ouput
            sc_dir = os.path.join(output_folder, scan_id)
            reg_folder = XnatUtils.makedir(sc_dir, subdir=False)
            ala_dir = os.path.join(reg_folder, 'REG_ALA')
            aff_dir = os.path.join(reg_folder, 'AFF')
            reg_dir = os.path.join(reg_folder, 'REG_F3D')
            cpp_dir = os.path.join(reg_folder, 'CPP')
            ala_folder = XnatUtils.makedir(ala_dir, subdir=False)
            aff_folder = XnatUtils.makedir(aff_dir, subdir=False)
            f3d_folder = XnatUtils.makedir(reg_dir, subdir=False)
            cpp_folder = XnatUtils.makedir(cpp_dir, subdir=False)
            self.time_writer("reg_aladin with ref %s and flo %s" %
                             (self.target['nii'], res_dict['nii']))
            aladin_output = os.path.join(ala_folder, "%s_2_%s_reg_aladin.nii" %
                                                     (scan_id, ARGS.target_id))
            affine_fpath = os.path.join(aff_folder,
                                        "%s_2_%s_affine_transformation.txt" %
                                        (scan_id, ARGS.target_id))

            cmd = REG_ALADIN_CMD.format(exe_path=self.reg_aladin_exe,
                                        ref=self.target['nii'],
                                        flo=res_dict['nii'],
                                        res=aladin_output,
                                        aff=affine_fpath,
                                        args=ARGS.args_reg_aladin)
            self.run_system_cmd(cmd)
            # Check that the affine file exists:
            if not os.path.exists(affine_fpath):
                err = 'Reg_aladin failed. File %s not found.'
                raise Exception(err % affine_fpath)

            # REG_F3D
            self.time_writer("reg_f3d with ref %s and flo %s and aff %s" %
                             (self.target['nii'],
                              res_dict['nii'],
                              affine_fpath))
            f3d_output = os.path.join(f3d_folder, "%s_2_%s_reg_f3d.nii" %
                                                  (scan_id, ARGS.target_id))
            f3d_cpp = os.path.join(cpp_folder, "%s_2_%s_reg_f3d_cpp.nii" %
                                               (scan_id, ARGS.target_id))
            cmd = REG_F3D_CMD.format(exe_path=self.reg_f3d_exe,
                                     ref=self.target['nii'],
                                     flo=res_dict['nii'],
                                     res=f3d_output,
                                     cpp=f3d_cpp,
                                     aff=affine_fpath,
                                     args=ARGS.args_regf3d)
            self.run_system_cmd(cmd)
            XnatUtils.gzip_nii(ala_folder)
            XnatUtils.gzip_nii(f3d_folder)
            XnatUtils.gzip_nii(cpp_folder)
            self.outputs[scan_id] = [{'label': 'reg_aladin_results',
                                      'image': aladin_output+'.gz'},
                                     {'label': 'reg_f3d_results',
                                      'image': f3d_output+'.gz'}]

            # Generate DICOM version of the reg_f3d results:
            convert_nifti_2_dicoms(f3d_output+'.gz', dcm_obj_sorted_list,
                                   self.sources[scan_id]['dcm'],
                                   os.path.join(output_folder, 'OSIRIX'),
                                   label=("%s_reg_f3d" % scan_id))

        # Make PDF
        self.make_pdf()

        # Zip the DICOMs output:
        initdir = os.getcwd()
        # Zip all the files in the directory
        zip_name = os.path.join(self.jobdir, 'outputs', 'OSIRIX', 'osirix.zip')
        os.chdir(os.path.join(self.jobdir, 'outputs', 'OSIRIX'))
        os.system('zip -r %s * > /dev/null' % zip_name)
        # return to the initial directory:
        os.chdir(initdir)

    def finish(self):
        """Method to copy the results in dax.RESULTS_DIR."""
        out_dir = os.path.join(self.jobdir, 'outputs')
        # Organise the outputs:
        ala_dir = os.path.join(out_dir, 'REG_ALA')
        aff_dir = os.path.join(out_dir, 'AFF')
        reg_dir = os.path.join(out_dir, 'REG_F3D')
        cpp_dir = os.path.join(out_dir, 'CPP')
        # Copy files:
        for scan_id, res_dict in self.sources.items():
            for folder in ['REG_ALA', 'REG_F3D', 'AFF', 'CPP']:
                old_path = glob.glob(os.path.join(out_dir, scan_id,
                                                  folder, '*'))
                new_path = os.path.join(out_dir, folder,
                                        os.path.basename(old_path))
                shutil.copy(old_path, new_path)
        # Zipping all the dicoms in the OSIRIX folder and keep the zip
        zip_osirix = os.path.join(out_dir, 'OSIRIX', 'osirix.zip')
        results_dict = {'PDF': self.pdf_final,
                        'REG_ALA': ala_dir,
                        'AFF': aff_dir,
                        'REG_F3D': reg_dir,
                        'CPP': cpp_dir,
                        'OSIRIX': zip_osirix}
        # Upload data:
        self.upload_dict(results_dict)
        self.end()

    def plot_images(self, fig, data, subplot_index, label):
        """Method to plot the images on the PDF for the first page.

        :param fig: figure from matplotlib
        :param data: numpy array of the images (3D) to display
        :param subplot_index: index of the line to display (4)
        :param label: y-label
        :return: None
        """
        ax = fig.add_subplot(4, 3, 3*subplot_index+1)
        ax.imshow(np.rot90(data[:, :, data.shape[2]/2]), cmap='gray')
        if subplot_index == 0:
            ax.set_title('Axial', fontsize=7)
        ax.set_ylabel(label, fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
        ax = fig.add_subplot(4, 3, 3*subplot_index+2)
        ax.imshow(np.rot90(data[:, data.shape[1]/2, :]), cmap='gray')
        if subplot_index == 0:
            ax.set_title('Coronal', fontsize=7)
        ax.set_axis_off()
        ax = fig.add_subplot(4, 3, 3*subplot_index+3)
        ax.imshow(np.rot90(data[data.shape[0]/2, :, :]), cmap='gray')
        if subplot_index == 0:
            ax.set_title('Sagittal', fontsize=7)
        ax.set_axis_off()

    @staticmethod
    def save_pdf_page(fig, pdf_path, date, pdf_page_number, pdf_pages):
        """Method to save the pdf page to a file.

        :param fig: figure from matplotlib
        :param pdf_path: numpy array of the images (3D) to display
        :param number_images: number of images for the subplot (6)
        :param subplot_index: index of the line to display (1-6)
        :param label: y-label
        :return: None
        """
        # Set footer and title
        fig.tight_layout()
        plt.figtext(0.5, 0.985,
                    '-- Registration Prostate Pipeline PDF report --',
                    horizontalalignment='center', fontsize=10)
        footer = 'Date: %s -- page %s/%s -- PDF generated by TIG laboratory \
at UCL, London' % (str(date), str(pdf_page_number), str(pdf_pages))
        plt.figtext(0.5, 0.02, footer, horizontalalignment='center',
                    fontsize=8)
        # Writing Figure
        fig.savefig(pdf_path, transparent=True, orientation='portrait',
                    dpi=100)
        plt.close(fig)

    def make_pdf(self):
        """Method to make the PDF for the spider.

        :return: None
        """
        # Define output files
        # Variables:
        date = datetime.datetime.now()
        page_count = 1

        # PDF path:
        pdf_pages_list = list()
        total_nb_pages = len(self.sources)
        ref_data = open_nifti(self.target['nii'])
        for scan_id, res_dict in self.sources.items():
            fig = plt.figure(page_count, figsize=(7.5, 10))
            self.plot_images(fig, ref_data, 0,
                             'Target - %s - %s' %
                             (self.target['ID'],
                              self.target['type']))
            source_data = open_nifti(res_dict['nii'])
            self.plot_images(fig, source_data, 1,
                             'Source - %s - %s' %
                             (self.sources[scan_id]['ID'],
                              self.sources[scan_id]['type']))

            for index, out_dict in enumerate(self.outputs[scan_id]):
                if not os.path.exists(out_dict['image']):
                    err = '%s output image not found.'
                    raise Exception(err % (out_dict['image']))

                # Open niftis with nibabel
                out_data = open_nifti(out_dict['image'])
                self.plot_images(fig, out_data, index+2, out_dict['label'])

            # Save figure to PDF page
            pdf_name = 'Reg_%s_2_%s_page_%s.pdf' % (ARGS.target_id, scan_id,
                                                    page_count)
            pdf_page = os.path.join(self.jobdir, pdf_name)
            self.save_pdf_page(fig, pdf_page, date, page_count, total_nb_pages)
            pdf_pages_list.append(pdf_page)
            page_count += 1

        # Join the two pages for the PDF:
        cmd = 'gs -q -sPAPERSIZE=letter -dNOPAUSE -dBATCH \
-sDEVICE=pdfwrite -sOutputFile=%s' % (self.pdf_final)
        for page in pdf_pages_list:
            cmd = '%s %s' % (cmd, page)
        self.time_writer('INFO:saving final PDF: %s ' % cmd)
        os.system(cmd)


def open_nifti(nifti_path):
    """Open the nifti from the path given to display for the PDF.

    :param nifti_path: path for the nifti
    :return: numpy array of data (3D matrix)
    """
    # Open niftis with nibabel
    f_img = nib.load(nifti_path)
    f_img_data = f_img.get_data()
    # Draw
    if len(f_img_data.shape) == 3:
        data = f_img_data
    elif len(f_img_data.shape) == 4:
        data = f_img_data[:, :, :, f_img_data.shape[3]/2]
    return data


def write_dicom(pixel_array, filename, ds_copy, ds_ori, volume_number,
                series_number, sop_id):
    """Write a dicom from a pixel_array (numpy).

    :param pixel_array: 2D numpy ndarray.
                        If pixel_array is larger than 2D, errors.
    :param filename: string name for the output file.
    :param ds_copy: pydicom object with the header that need to be copy
    :param ds_ori: original pydicom object of the pixel_array
    :param volume_number: number of the volume being processed
    :param series_number: number of the series being processed
    :param sop_id: SOPInstanceUID for the DICOM
    """
    # Set to zero negatives values in the image:
    pixel_array[pixel_array < 0] = 0

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
    ds.SeriesDescription = ds_ori.SeriesDescription + ' reg_f3d'
    sop_uid = sop_id + str(datetime.datetime.now()).replace('-', '')\
                                                   .replace(':', '')\
                                                   .replace('.', '')\
                                                   .replace(' ', '')
    ds.SOPInstanceUID = sop_uid[:-1]
    ds.ProtocolName = ds_ori.ProtocolName
    ds.InstanceNumber = volume_number+1

    # Copy from T2 the orientation tags:
    for tag in TAGS_TO_COPY:
        if tag in ds_copy:
            ds[tag] = ds_copy[tag]

    # Set the Image pixel array
    if pixel_array.dtype != np.uint16:
        pixel_array = pixel_array.astype(np.uint16)
    ds.PixelData = pixel_array.tostring()

    # Save the image
    ds.save_as(filename)


def convert_nifti_2_dicoms(nifti_path, dcm_targets, dicom_source,
                           output_folder, label=None):
    """Convert 4D niftis generated by reg_f3d into DICOM files.

    :param nifti_path: path to the nifti file
    :param dcm_targets: list of pydicom object from the target image for
                        header info
    :param dicom_source: one dicom file from the source for header info
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
    if not os.path.isfile(dicom_source):
        err = "DICOM File %s not found after reg_f3d."
        raise Exception(err % dicom_source)
    adc_dcm_obj = dicom.read_file(dicom_source,
                                  force=True)

    # Make output_folder:
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Series Number and SOP UID
    ti = time.time()
    series_number = 86532 + int(str(ti)[2:4]) + int(str(ti)[4:6])
    sop_id = adc_dcm_obj.SOPInstanceUID.split('.')
    sop_id = '.'.join(sop_id[:-1])+'.'

    for volume_index in range(f_img_data.shape[2]):
        if f_img_data.shape[2] > 100:
            filename = os.path.join(output_folder, '%s_%03d.dcm' %
                                                   (label, volume_index+1))
        elif f_img_data.shape[2] > 10:
            filename = os.path.join(output_folder, '%s_%02d.dcm' %
                                                   (label, volume_index+1))
        else:
            filename = os.path.join(output_folder, '%s_%d.dcm' %
                                                   (label, volume_index+1))
        write_dicom(np.rot90(f_img_data[:, :, volume_index]), filename,
                    dcm_targets[volume_index], adc_dcm_obj, volume_index,
                    series_number, sop_id)

if __name__ == '__main__':
    ARGS = parse_args()
    # generate spider object:
    spider_obj = Spider_Registration_Prostate(spider_path=sys.argv[0],
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
