"""Spider to register VERDICT modalities for prostate cancer project.

Author:         Benjamin Yvernault
contact:        b.yvernault@ucl.ac.uk
Spider name:    Registration_Prostate
Spider version: 1.0.0
Creation date:  2016-03-15 13:33:03.911277
Purpose:        Register VERDICT modalities
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
from dax import spiders, XnatUtils, SessionSpider

__author__ = "Benjamin Yvernault"
__email__ = "b.yvernault@ucl.ac.uk"
__purpose__ = "Register VERDICT modalities."
__spider_name__ = "Registration_Verdict"
__version__ = "1.0.0"
__modifications__ = """2016-06-21 16:41:03 - Original write
"""

REG_ALADIN_CMD = "{exe_path} \
-ref {ref} \
-flo {flo} \
-res {res} \
-aff {aff} \
{args}"
REG_RESAMPLE_CMD = "{exe_path} \
-ref {ref} \
-flo {flo} \
-aff {aff} \
-res {res} \
{args}"
DEFAULT_ARGS_REG_ALADIN = " -rigOnly "
DEFAULT_ARGS_REG_RESAMPLE = ""
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
        --regResample : path to reg_resample's executable
        --argRegResample : arguments for reg_resample.
          Default: -ln 4 -lp 4 -jl 0.1 -be 0.05 -maxit 250 -lncc 0 5.0 -sx 2.5
        --openmp_core : number of core use by reg_aladin. Default: one

    :return: argument parser object created by parse_args()
    """
    ap = spiders.get_session_argparser("Spider_%s" % __spider_name__,
                                       __purpose__)
    ap.add_argument("--scansid", dest="scans_id", required=True,
                    help="VERDICT Scans ID from XNAT.")
    ap.add_argument("--regAladin", dest="reg_aladin_exe", required=True,
                    help="path to reg_aladin's executable.")
    ap.add_argument("--argsRegAladin", dest="args_reg_aladin",
                    help="Argument for reg_aladin. \
Default: -rigOnly.", default=DEFAULT_ARGS_REG_ALADIN)
    ap.add_argument("--regResample", dest="reg_resample_exe", required=True,
                    help="path to reg_resample's executable.")
    ap.add_argument("--argRegResample", dest="args_reg_resample",
                    help="Argument for reg_resample. \
Default: None.", default=DEFAULT_ARGS_REG_RESAMPLE)
    ap.add_argument("--openmp_core", dest="openmp_core", default=1,
                    help="Number of core used by reg_aladin.")
    return ap.parse_args()


class Spider_Registration_Verdict(SessionSpider):
    """Session Spider class to do: Register VERDICT modalities.

    :param spider_path: spider file path
    :param jobdir: directory for temporary files
    :param scans_id: VERDICT scans ID XNAT
    :param xnat_project: project ID on XNAT
    :param xnat_subject: subject label on XNAT
    :param xnat_session: experiment label on XNAT
    :param reg_aladin_exe: reg_aladin executable path
    :param reg_resample_exe: reg_resample executable path
    :param xnat_host: host for XNAT if not set in environment variables
    :param xnat_user: user for XNAT if not set in environment variables
    :param xnat_pass: password for XNAT if not set in environment variables
    :param suffix: suffix to the assessor creation
    """

    def __init__(self, spider_path, jobdir, scans_id,
                 xnat_project, xnat_subject, xnat_session,
                 reg_aladin_exe='reg_aladin', reg_resample_exe='reg_resample',
                 args_reg_aladin=DEFAULT_ARGS_REG_ALADIN,
                 args_reg_resample=DEFAULT_ARGS_REG_RESAMPLE,
                 xnat_host=None, xnat_user=None, xnat_pass=None, suffix=""):
        """Entry point for Spider_Registration_Prostate Class."""
        super(Spider_Registration_Verdict, self).__init__(spider_path, jobdir,
                                                          xnat_project,
                                                          xnat_subject,
                                                          xnat_session,
                                                          xnat_host,
                                                          xnat_user,
                                                          xnat_pass, suffix)
        # Inputs
        self.acquisitions = dict()
        self.dicom = ''
        self.scans_id = XnatUtils.get_input_list(scans_id, list())

        # Outputs
        self.pdf_final = os.path.join(self.jobdir, 'Registration_VERDICT.pdf')

        # Check Executable:
        self.reg_aladin_exe = self.check_executable(reg_aladin_exe,
                                                    'reg_aladin')
        self.reg_resample_exe = self.check_executable(reg_resample_exe,
                                                      'reg_resample')
        self.args_reg_aladin = args_reg_aladin
        self.args_reg_resample = args_reg_resample

    def pre_run(self):
        """Method to download data from XNAT.

        :param argument_parse: argument parser object return by parse_args()
        """
        # Make directory
        input_folder = XnatUtils.makedir(os.path.join(self.jobdir, 'inputs'),
                                         subdir=False)

        # Download scans:
        self.time_writer('Connection to XNAT')
        xnat = XnatUtils.get_interface(host=self.host,
                                       user=self.user,
                                       pwd=self.pwd)
        # Download one DICOM
        dcm_dir = XnatUtils.makedir(os.path.join(input_folder, 'DICOM'),
                                    subdir=False)
        scan = XnatUtils.select_obj(xnat,
                                    self.xnat_project,
                                    self.xnat_subject,
                                    self.xnat_session,
                                    self.scans_id[0])
        sdcm_obj = scan.resource('DICOM')
        self.time_writer('Downloading DICOM for scan ID %s ...'
                         % self.scans_id[0])
        self.dicom = XnatUtils.download_file_from_obj(dcm_dir, sdcm_obj)

        # Download NIFTIs
        index = 1
        for scan_id in self.scans_id:
            scan_info = {}
            scan_dir = XnatUtils.makedir(os.path.join(input_folder, scan_id),
                                         subdir=False)

            self.time_writer('Downloading scan ID %s ...' % scan_id)
            scan = XnatUtils.select_obj(xnat,
                                        self.xnat_project,
                                        self.xnat_subject,
                                        self.xnat_session,
                                        scan_id)
            snii_obj = scan.resource('NIFTI')
            scan_info['4D'] = XnatUtils.download_file_from_obj(scan_dir,
                                                               snii_obj)
            scan_info['type'] = scan.attrs.get('type')
            scan_info['ID'] = scan_id
            if 'b3000' in scan_info['type'].lower() and \
               len([o for o in self.acquisitions.get(index, list())
                    if 'b3000' in o['type'].lower()]) > 0:
                index += 1

            if index in self.acquisitions:
                self.acquisitions[index].append(scan_info)
            else:
                self.acquisitions[index] = [scan_info]
        xnat.disconnect()
        self.time_writer('Disconnection of XNAT')

    def run(self):
        """Method running the process for the spider on the inputs data."""
        output_folder = XnatUtils.makedir(os.path.join(self.jobdir, 'outputs'))
        dcm_folder = XnatUtils.makedir(os.path.join(output_folder, 'DICOM'))

        for i in range(1, len(self.acquisitions.keys()) + 1):
            for index, scan_info in enumerate(self.acquisitions[i]):
                # Step 1:
                # Transform the 4D Nifti into 3D images for registration:
                self.time_writer('Splitting nifti %s ...' % scan_info['ID'])
                self.acquisitions[i][index]['3D'] = split_nifti_4D_3Ds(
                                                        scan_info['4D'])

                if 'b3000' not in scan_info['type'].lower():
                    outdir = XnatUtils.makedir(os.path.join(output_folder,
                                                            scan_info['ID']))
                    # Step 2-3-4 see register nifti
                    self.time_writer('Registration nifti ...')

                    nii_reg = self.register_nifti(
                                        self.acquisitions[i][index],
                                        self.acquisitions[i][index-1],
                                        outdir, i)

                    self.acquisitions[i][index]['reg'] = nii_reg

                    # Generate DICOM version of the reg_f3d results:
                    convert_nifti_2_dicoms(nii_reg,
                                           self.dicom,
                                           self.dicom,
                                           dcm_folder,
                                           scan_info['type'],
                                           label=("%s_%s_reg"
                                                  % (scan_info['ID'],
                                                     scan_info['type'])))
        shutil.move(self.dicom, dcm_folder)

        # Make PDF
        # self.make_pdf()

        # Zip the DICOMs output:
        """initdir = os.getcwd()
        # Zip all the files in the directory
        zip_name = os.path.join(self.jobdir, 'outputs', 'OSIRIX', 'osirix.zip')
        os.chdir(os.path.join(self.jobdir, 'outputs', 'OSIRIX'))
        os.system('zip -r %s * > /dev/null' % zip_name)
        # return to the initial directory:
        os.chdir(initdir)"""

    def finish(self):
        """Method to copy the results in dax.RESULTS_DIR."""
        """out_dir = os.path.join(self.jobdir, 'outputs')
        # Organise the outputs:
        ala_dir = XnatUtils.makedir(os.path.join(out_dir, 'REG_ALA'))
        aff_dir = XnatUtils.makedir(os.path.join(out_dir, 'AFF'))
        reg_dir = XnatUtils.makedir(os.path.join(out_dir, 'REG_F3D'))
        cpp_dir = XnatUtils.makedir(os.path.join(out_dir, 'CPP'))
        # Copy files:
        for scan_id, res_dict in self.sources.items():
            for folder in ['REG_ALA', 'REG_F3D', 'AFF', 'CPP']:
                old_path = glob.glob(os.path.join(out_dir, scan_id,
                                                  folder, '*'))[0]
                new_path = os.path.join(out_dir, folder)
                shutil.copy(old_path, new_path)
        # Zipping all the dicoms in the OSIRIX folder and keep the zip
        zip_osirix = os.path.join(out_dir, 'OSIRIX', 'osirix.zip')
        results_dict = {'PDF': self.pdf_final,
                        'ACQ1': ala_dir,
                        'ACQ2': aff_dir,
                        'OSIRIX': zip_osirix}
        # Upload data:
        self.upload_dict(results_dict)
        self.end()"""

    def make_pdf(self):
        """Method to make the PDF for the spider.

        :return: None
        """
        print 'PDF'

    def register_nifti(self, source_info, target_info, output_folder,
                       acquisition_number):
        """Register the nifti source to the target.

        :param source_info: dictionary information of source image
        :param target_info: dictionary information of target image
        :param output_folder: path to the output folder
        :param acquisition_number: index of the acquisition
        :return: path to the 4D nifti register
        """
        # Variables:
        ala_dir = XnatUtils.makedir(os.path.join(output_folder, 'REG_ALA'))
        reg_dir = XnatUtils.makedir(os.path.join(output_folder, 'REG_RES'))
        b1tob0 = os.path.join(ala_dir, 'b1tob0.txt')
        b0_nii = os.path.join(ala_dir, 'b0_reg.nii')
        volume_niis = {0: b0_nii}
        # Step 2:
        # Register each scan b0 to the previous one
        # (e.g: b3000 <- b2000)
        cmd = REG_ALADIN_CMD.format(exe_path=self.reg_aladin_exe,
                                    ref=target_info['3D'][0],
                                    flo=source_info['3D'][0],
                                    res=b0_nii,
                                    aff=b1tob0,
                                    args=self.args_reg_aladin)
        self.run_system_cmd(cmd)

        # Step 3:
        # Propagate the transformation to the rest of the volume:
        for index, volume in enumerate(source_info['3D']):
            if index != 0:
                vol_nii = os.path.join(reg_dir, 'volume_%s_reg.nii' % index)
                cmd = REG_ALADIN_CMD.format(exe_path=self.reg_resample_exe,
                                            ref=target_info['3D'][0],
                                            flo=volume,
                                            res=vol_nii,
                                            aff=b1tob0,
                                            args=self.args_reg_resample)
                self.run_system_cmd(cmd)
                volume_niis[index] = vol_nii

        # Step 4:
        # Put back the nifti together as one volume
        final_nii = os.path.join(output_folder, '%s_%s_%d_reg.nii'
                                                % (source_info['ID'],
                                                   source_info['type'],
                                                   acquisition_number))
        join_nifti_3Ds_4D(volume_niis, final_nii)
        return final_nii

    def generate_big_nifti(self, nifti_path):
        """Generate big nifti with all VERDICT acquisition files.

        :param nifti_path: nifti path
        """
        for i in range(1, len(self.acquisitions.keys()) + 1):
            f_img = nib.load(self.acquisitions[i][0]['4D'])
            f_img_data = f_img.get_data()
            data = np.zeros(shape=(f_img_data.shape[0],
                                   f_img_data.shape[1],
                                   f_img_data.shape[2],
                                   f_img_data.shape[3],
                                   len(self.acquisitions[i])))
            for index, scan_info in enumerate(self.acquisitions[i]):
                if index == 0:
                    key = '4D'
                else:
                    key = 'reg'

                # Open niftis with nibabel
                f_img = nib.load(scan_info[key])
                f_img_data = f_img.get_data()
                # Draw
                data[:, :, :, :, index] = f_img_data
                nii_5d = nib.Nifti1Image(data, affine=f_img.affine)
            acq_dir = XnatUtils.makedir(os.path.join(self.jobdir,
                                                     'outputs',
                                                     'ACQ%d' % index))
            filename = '%s_acquisition%d.nii' % (self.xnat_session, index)
            nii_file = os.path.join(acq_dir, filename)
            nib.save(nii_5d, nii_file)


def split_nifti_4D_3Ds(nifti_path):
    """Split 4D niftis into 3Ds.

    :param nifti_path: path for the nifti
    :return: list of path of new 3D niftis
    """
    # list of images:
    li_niis = list()
    # Open niftis with nibabel
    f_img = nib.load(nifti_path)
    f_img_data = f_img.get_data()
    # Draw
    if len(f_img_data.shape) < 3:
        return [nifti_path]
    elif len(f_img_data.shape) > 4:
        print 'Warning: images dimension > 4. Can not be split.'
    else:
        for index in range(f_img_data.shape[3]):
            data = f_img_data[:, :, :, index]
            nii_3d = nib.Nifti1Image(data, affine=f_img.affine)
            nii_file = '%s_%d.nii' % (os.path.splitext(nifti_path)[0],
                                      index)
            nib.save(nii_3d, nii_file)
            li_niis.append(nii_file)
    return li_niis


def join_nifti_3Ds_4D(li_nii, nifti_path):
    """Join 3D niftis into 4D.

    :param li_nii: dictionary of nifti file paths in order
    :param nifti_path: nifti_path
    """
    f_img = nib.load(li_nii[0])
    f_img_data = f_img.get_data()
    data = np.zeros(shape=(f_img_data.shape[0],
                           f_img_data.shape[1],
                           f_img_data.shape[2],
                           len(li_nii.keys())))
    for index, nii in li_nii.items():
        # Open niftis with nibabel
        f_img = nib.load(nii)
        f_img_data = f_img.get_data()
        # Draw
        if len(f_img_data.shape) != 3:
            print 'Warning: nifti image not 3D. Can not be join.'
        else:
            data[:, :, :, index] = f_img_data
    nii_4d = nib.Nifti1Image(data, affine=f_img.affine)
    nib.save(nii_4d, nifti_path)


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
                series_number, sop_id, stype):
    """Write a dicom from a pixel_array (numpy).

    :param pixel_array: 2D numpy ndarray.
                        If pixel_array is larger than 2D, errors.
    :param filename: string name for the output file.
    :param ds_copy: pydicom object with the header that need to be copy
    :param ds_ori: original pydicom object of the pixel_array
    :param volume_number: number of the volume being processed
    :param series_number: number of the series being processed
    :param sop_id: SOPInstanceUID for the DICOM
    :param stype: type of the scan
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
    if 'SeriesDescription' not in ds:
        ds.SeriesDescription = stype + ' registered'
    else:
        ds.SeriesDescription = ds_ori.SeriesDescription + ' registered'
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
                           output_folder, stype, label=None):
    """Convert 4D niftis generated by reg_f3d into DICOM files.

    :param nifti_path: path to the nifti file
    :param dcm_targets: list of pydicom object from the target image for
                        header info
    :param dicom_source: one dicom file from the source for header info
    :param output_folder: folder where the DICOM files will be saved
    :param stype: type of the scan
    :param label: name for the output dicom files
    :return: None
    """
    if not os.path.isfile(nifti_path):
        raise Exception("File %s not found after registration." % nifti_path)
    # Load image from NIFTI
    f_img = nib.load(nifti_path)
    f_img_data = f_img.get_data()

    # Load dicom headers
    if not os.path.isfile(dicom_source):
        err = "DICOM File %s not found after registration."
        raise Exception(err % dicom_source)
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
                    series_number, sop_id, stype)

if __name__ == '__main__':
    ARGS = parse_args()
    # generate spider object:
    spider_obj = Spider_Registration_Verdict(
                        spider_path=sys.argv[0],
                        jobdir=ARGS.temp_dir,
                        scans_id=ARGS.scans_id,
                        xnat_project=ARGS.proj_label,
                        xnat_subject=ARGS.subj_label,
                        xnat_session=ARGS.sess_label,
                        reg_aladin_exe=ARGS.reg_aladin_exe,
                        reg_resample_exe=ARGS.reg_resample_exe,
                        args_reg_aladin=ARGS.args_reg_aladin,
                        args_reg_resample=ARGS.args_reg_resample,
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
    # spider_obj.finish()
