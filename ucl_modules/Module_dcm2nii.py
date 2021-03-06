"""Module to generate NIFTI from DICOM with dcm2nii."""

from dax import XnatUtils, ScanModule
import os
import glob
import logging
import nibabel as nib
import shutil
import subprocess as sb

LOGGER = logging.getLogger('dax')

DEFAULT_TPM_PATH = '/tmp/dcm2nii_temp/'
DEFAULT_MODULE_NAME = 'dcm2nii'
DEFAULT_TEXT_REPORT = 'ERROR/WARNING for dcm2nii :\n'
DCM2NII_CMD = '''{dcm2nii} -a n -e n -d n -g y -f n -n y -p n \
-v y -x n -r n {dicom}'''
AVOID_SCANTYPES = ['PhoenixZIPReport']
DCMDJPEG_TEMPLATE = """{dcmdjpeg} {original_dcm} {new_dcm} > /dev/null"""


class Module_dcm2nii(ScanModule):
    """Module to convert dicom to nifti using dcm2nii."""

    def __init__(self, mod_name=DEFAULT_MODULE_NAME,
                 directory=DEFAULT_TPM_PATH, email=None,
                 text_report=DEFAULT_TEXT_REPORT,
                 zip_dicoms=False,
                 dcm2nii_exe='dcm2nii',
                 dcmdjpeg_exe='dcmdjpeg'):
        """init function overridden from base-class."""
        super(Module_dcm2nii, self).__init__(mod_name, directory, email,
                                             text_report=text_report)
        self.dcm2nii_exe = check_executable(dcm2nii_exe)
        self.dcmdjpeg_exe = check_executable(dcmdjpeg_exe)
        self.dicom_paths = list()
        self.zip_dicoms = zip_dicoms

    def prerun(self, settings_filename=''):
        """prerun function overridden from base-class."""
        # make directory
        self.make_dir(settings_filename)

    def afterrun(self, xnat, project):
        """afterrun function overridden from base-class."""
        # send report
        if self.send_an_email:
            self.send_report()

        # clean the directory created
        try:
            os.rmdir(self.directory)
        except Exception:
            LOGGER.warn('dcm2nii -- afterrun -- %s not empty. Could not \
delete it.' % self.directory)

    def needs_run(self, cscan, xnat):
        """needs_run function overridden from base-class.

        cscan = CacheScan object from XnatUtils
        return True or False
        """
        # Check output
        if XnatUtils.has_resource(cscan, 'NIFTI'):
            LOGGER.debug('Has NIFTI')
            return False
        # Check input
        if not XnatUtils.has_resource(cscan, 'DICOM'):
            LOGGER.debug('no DICOM resource')
            return False

        return True

    def run(self, scan_info, scan_obj):
        """run function to convert dicom to parrec to nifti and upload data."""
        # clean tmp folder
        XnatUtils.clean_directory(self.directory)

        if not len(scan_obj.resource('DICOM').files().get()) > 0:
            LOGGER.debug('no DICOM files')
        elif scan_info['type'] in AVOID_SCANTYPES:
            LOGGER.info('avoid this scan type: {}'.format(scan_info['type']))
        else:
            LOGGER.debug('downloading all DICOMs...')

            self.dicom_paths = XnatUtils.download_files_from_obj(
                self.directory, scan_obj.resource('DICOM'))

            if not self.dicom_paths:
                msg = """dcm2nii -- %s -- No proper DICOM found in \
    resource DICOM on XNAT"""
                LOGGER.error(msg % scan_info['scan_id'])
                msg = 'No proper DICOM found in resource DICOM on XNAT'
                self.log_warning_error(msg, scan_info, error=True)
            else:
                # convert dcm to nii
                dcm_dir = os.path.dirname(self.dicom_paths[0])

                # ZIP the DICOM if more than one
                if len(self.dicom_paths) > 1 and self.zip_dicoms:
                    self.zipping_dicoms(scan_obj, dcm_dir)

                # if only one DICOM and it's a zip, unzip
                if len(self.dicom_paths) == 1 and \
                   self.dicom_paths[0].endswith('.zip'):
                    dcm_dir = os.path.dirname(self.dicom_paths[0])
                    os.system('unzip -d %s -j %s > /dev/null'
                              % (dcm_dir, self.dicom_paths[0]))
                    os.remove(self.dicom_paths[0])
                    self.dicom_paths = self.get_dicom_list(dcm_dir)

                if not self.dcm2nii(self.dicom_paths[0]):
                    # Convert dcm via dcmdjpeg
                    dicom_paths_djpeg = self.dcmdjpeg()
                    # try again dcm2nii
                    self.dcm2nii(dicom_paths_djpeg[0])
                    dcm_dir = os.path.dirname(dicom_paths_djpeg[0])

                # Check if Nifti created:
                nifti_list = [
                    f for f in os.listdir(dcm_dir)
                    if f.endswith('.nii.gz') or f.endswith('.nii')]
                if not nifti_list:
                    msg = "dcm2nii -- %s -- DCM --> NII ( preprocess \
dicom with dcmdjpeg ) conversion failure"
                    LOGGER.warn(msg % scan_info['scan_id'])
                    msg = 'Fail to convert DICOM to NIFTI '
                    self.log_warning_error(msg, scan_info)
                else:
                    # UPLOADING THE RESULTS
                    self.upload_converted_images(
                        dcm_dir, scan_obj, scan_info)

            # clean tmp folder
            LOGGER.debug('clean temp directory...')
            XnatUtils.clean_directory(self.directory)

    @staticmethod
    def is_dicom(fpath):
        """check if the file is a DICOM medical data.

        :param fpath: path of the file
        :return boolean: true if it's a DICOM, false otherwise
        """
        file_call = '''file {fpath}'''.format(fpath=fpath)
        output = sb.check_output(file_call.split())
        if 'dicom' in output.lower():
            return True

        return False

    def get_dicom_list(self, directory):
        """get the list of DICOMs file from the directory.

        :param directory: directory containing the DICOM files.
        :return list(): list of filepaths that are dicoms in directory
        """
        fnames = os.listdir(directory)
        dicom_paths = list()
        for fname in fnames:
            fpath = os.path.join(directory, fname)
            if self.is_dicom(fpath):
                dicom_paths.append(fpath)

        return dicom_paths

    @staticmethod
    def zipping_dicoms(scan_obj, dcm_dir):
        """Zipping the dicoms.

        :param scan_obj: scan pyxnat object
        :param dcm_dir: folder containing the dicoms
        """
        LOGGER.debug('   --> more than one \
dicom files, zipping dicoms.')
        fzip = 'dicoms.zip'
        initdir = os.getcwd()
        # Zip all the files in the directory
        os.chdir(dcm_dir)
        os.system('zip -r {} * > /dev/null'.format(fzip))
        fzip_path = os.path.join(dcm_dir, fzip)
        # return to the initial directory:
        os.chdir(initdir)
        # upload
        if os.path.exists(fzip_path):
            LOGGER.debug('   --> uploading zip dicoms')
            scan_obj.resource('DICOM').delete()
            scan_obj.resource('DICOM').put_zip(fzip_path,
                                               overwrite=True,
                                               extract=False)

    def dcm2nii(self, dicom_path):
        """convert dicom to nifti using dcm2nii."""
        LOGGER.debug('convert dcm to nii...')
        dcm2nii_cmd = DCM2NII_CMD.format(dcm2nii=self.dcm2nii_exe,
                                         dicom=dicom_path)
        try:
            sb.check_output(dcm2nii_cmd.split())
        except sb.CalledProcessError:
            LOGGER.debug("DCM --> NII conversion failed")
            return False

        return True

    def dcmdjpeg(self):
        """Decompress the dicom from jpeg.

        :return: list of dicoms
        """
        LOGGER.debug('run dcmdjpeg on the DICOMs.')
        dicom_paths = []
        dcm_dir = os.path.join(os.path.dirname(self.dicom_paths[0]),
                               'DCMDJPEGEDs')
        if not os.path.exists(dcm_dir):
            os.makedirs(dcm_dir)
        for dicom in self.dicom_paths:
            root, ext = os.path.splitext(dicom)
            dcm_p = os.path.join(dcm_dir, os.path.basename(root))
            new_dicom = "%s_dcmdjpeged.%s" % (dcm_p, ext)
            dcmdjpeg_cmd = DCMDJPEG_TEMPLATE.format(
                dcmdjpeg=self.dcmdjpeg_exe,
                original_dcm=dicom,
                new_dcm=new_dicom)
            os.system(dcmdjpeg_cmd)
            dicom_paths.append(new_dicom)
        return dicom_paths

    def upload_converted_images(self, dcm_dir, scan_obj, scan_info):
        """upload the images after checking them."""
        # Local variables
        nifti_list = []
        bval_fpath = ''
        bvec_fpath = ''

        # G et the bvec/bval files and NIFTI from the folder:
        for fpath in glob.glob(os.path.join(dcm_dir, '*')):
            if os.path.isfile(fpath):
                if fpath.lower().endswith('.bval'):
                    bval_fpath = fpath
                if fpath.lower().endswith('.bvec'):
                    bvec_fpath = fpath
                if fpath.lower().endswith('.nii.gz'):
                    nifti_list.append(fpath)
                if fpath.lower().endswith('.nii'):
                    os.system('gzip ' + fpath)
                    nifti_list.append(fpath + '.gz')

        # Check NIFTI:
        good_to_upload = self.check_outputs(scan_info, nifti_list,
                                            bval_fpath, bvec_fpath)
        # Upload files:
        if good_to_upload:
            if os.path.isfile(bval_fpath) and os.path.isfile(bvec_fpath):
                # BVAL/BVEC
                XnatUtils.upload_file_to_obj(bval_fpath,
                                             scan_obj.resource('BVAL'),
                                             remove=True)
                XnatUtils.upload_file_to_obj(bvec_fpath,
                                             scan_obj.resource('BVEC'),
                                             remove=True)
                # keep the NII with the same name than the BVAL/BVEC
                nifti_list = filter(lambda x: x[:-7] == bval_fpath[:-5],
                                    nifti_list)
                XnatUtils.upload_files_to_obj(nifti_list,
                                              scan_obj.resource('NIFTI'),
                                              remove=True)
            else:
                # NII
                XnatUtils.upload_files_to_obj(nifti_list,
                                              scan_obj.resource('NIFTI'),
                                              remove=True)

            # more than one NIFTI uploaded
            if len(nifti_list) > 1:
                LOGGER.warn('''dcm2nii -- %s -- more than one NIFTI upload'''
                            % scan_info['scan_id'])
                self.log_warning_error('more than one NIFTI upload', scan_info)

    def check_outputs(self, scan_info, nifti_list, bval, bvec):
        """Check that the outputs are right (opening nifti works).

        :param scan_info: scan information dictionary
        :param nifti_list: python list of nifti paths
        :param bval: python list of bval paths
        :param bvec: python list of bvec paths
        :return boolean: true if outputs are fine, false otherwise
        """
        for nifti_fpath in nifti_list:
            try:
                nib.load(nifti_fpath)
            except Exception:
                LOGGER.warn("dcm2nii -- %s -- %s is not a proper NIFTI"
                            % (scan_info['scan_id'],
                               os.path.basename(nifti_fpath)))
                self.log_warning_error('check the scan DICOM. No upload. \
a non-valid nifti was created', scan_info, error=True)
                return False
        if bval or bvec:
            if not os.path.isfile(bval) or not os.path.isfile(bvec):
                LOGGER.warn("dcm2nii -- %s -- no bval/bvec generated \
(diffusion scan)" % scan_info['scan_id'])
                self.log_warning_error('check the scan DICOM. No upload. \
no bval/bvec were generated', scan_info, error=True)
                return False
        return True


def no_executable(cmd):
    try:
        return shutil.which(cmd) is None
    except AttributeError:  # Python < 3.3
        return not any(
            [os.path.isfile(os.path.join(path, cmd)) and
             os.access(os.path.join(path, cmd), os.X_OK)
             for path in os.environ["PATH"].split(os.pathsep)])


def check_executable(executable):
    """Method to check the executable.

    :param executable: executable path
    :param name: name of Executable
    :return: Complete path to the executable
    """
    if no_executable(executable):
        raise Exception("Executable '%s' not found" % (executable))

    pversion = sb.Popen([executable, '--version'],
                        stdout=sb.PIPE,
                        stderr=sb.PIPE)
    nve_version, _ = pversion.communicate()
    name = os.path.basename(executable)
    LOGGER.debug('%s version: %s' % (name, nve_version.strip().split('\n')[0]))
    return executable
