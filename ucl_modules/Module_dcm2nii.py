""" Module to generate NIFTI from DICOM with dcm2nii """
from dax import XnatUtils, ScanModule
from VUIIS_path_settings import DCM2NII_PATH, DCMDJPEG_PATH
import os
import glob
import dicom
import logging
import numpy as np
import nibabel as nib
import subprocess as sb

LOGGER = logging.getLogger('dax')

DEFAULT_TPM_PATH = '/tmp/dcm2nii_temp/'
DEFAULT_MODULE_NAME = 'dcm2nii'
DEFAULT_TEXT_REPORT = 'ERROR/WARNING for dcm2nii :\n'

class Module_dcm2nii(ScanModule):
    """ Module to convert dicom to nifti using dcm2nii """
    def __init__(self, mod_name=DEFAULT_MODULE_NAME, directory=DEFAULT_TPM_PATH, email=None, text_report=DEFAULT_TEXT_REPORT,
                 dcm2niipath=DCM2NII_PATH, dcmdjpegpath=DCMDJPEG_PATH):
        """ init function overridden from base-class"""
        super(Module_dcm2nii, self).__init__(mod_name, directory, email, text_report=text_report)
        self.dcm2niipath = dcm2niipath
        self.dcmdjpegpath = dcmdjpegpath
        self.dicom_paths = list()

    def prerun(self, settings_filename=''):
        """ prerun function overridden from base-class"""
        #make directory
        self.make_dir(settings_filename)

    def afterrun(self, xnat, project):
        """ afterrun function overridden from base-class"""
        #send report
        if self.send_an_email:
            self.send_report()

        #clean the directory created
        try:
            os.rmdir(self.directory)
        except:
            LOGGER.warn('dcm2nii -- afterrun -- '+self.directory+' not empty. Could not delete it.')

    def needs_run(self, cscan, xnat):
        """ needs_run function overridden from base-class
                cscan = CacheScan object from XnatUtils
            return True or False
        """
        #Variables:
        scan_info = cscan.info()

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
        """ run function to convert dicom to parrec to nifti and upload data"""
        if not len(scan_obj.resource('DICOM').files().get()) > 0:
            LOGGER.debug('no DICOM files')
        else:
            LOGGER.debug('downloading all DICOMs...')
            scan_obj.resource('DICOM').get(self.directory, extract=True)

            dcm_dir = os.path.join(self.directory, 'DICOM')
            self.set_dicom_list(dcm_dir)
            if self.dicom_paths:
                # Check for duplicate dicoms:
                self.check_duplicate_slices_dicom(scan_info)
                # convert dcm to nii
                conversion_status = self.dcm2nii(self.dicom_paths[0])

                if not conversion_status:
                    #Convert dcm via dcmdjpeg
                    self.dcmdjpeg(dcm_dir)
                    #try again dcm2nii
                    dcm_fpath = os.path.join(dcm_dir, 'final_1.dcm')
                    conversion_status = self.dcm2nii(dcm_fpath)

                #Check if Nifti created:
                nifti_list = [fname for fname in os.listdir(dcm_dir) if fname.endswith('.nii.gz') or fname.endswith('.nii')]
                if not nifti_list:
                    LOGGER.warn('''dcm2nii -- {scan} -- DCM --> NII ( preprocess dicom with dcmdjpeg ) conversion failure'''.format(scan=scan_info['scan_id']))
                    self.log_warning_error('Fail to convert DICOM to NIFTI ', scan_info)
                else:
                    ###### UPLOADING THE RESULTS ######
                    self.upload_converted_images(dcm_dir, scan_obj, scan_info)
            else:
                LOGGER.error('''dcm2nii -- {scan} -- No proper DICOM found in resource DICOM on XNAT'''.format(scan=scan_info['scan_id']))
                self.log_warning_error('No proper DICOM found in resource DICOM on XNAT', scan_info, error=True)

            # clean tmp folder
            self.clean_directory()

    @staticmethod
    def is_dicom(fpath):
        """
            check if the file is a DICOM medical data

            :param fpath: path of the file
            :return boolean: true if it's a DICOM, false otherwise
        """
        file_call = '''file {fpath}'''.format(fpath=fpath)
        output = sb.check_output(file_call.split())
        if 'dicom' in output.lower():
            return True

        return False

    def set_dicom_list(self, directory):
        """
            get the list of DICOMs file from the directory

            :param directory: directory containing the DICOM files.
        """
        fnames = os.listdir(directory)
        for fname in fnames:
            fpath = os.path.join(directory, fname)
            if self.is_dicom(fpath):
                self.dicom_paths.append(fpath)

    def check_duplicate_slices_dicom(self, scan_info):
        """
            Check for duplicate slices in the dicom

            :param dicom_fpaths: list of dicom files to convert
            :param scan_info: dictionary on scan from XNAT
            :return boolean: true if the dicoms are fine, false otherwise
        """
        if len(self.dicom_paths) > 2: #more than one dicom files downloaded in the folder
            new_size = [len(self.dicom_paths), 3]
            orien_mat = np.zeros(new_size)

            for index, dicom_path in enumerate(self.dicom_paths):
                #read the DICOM header
                dcm_header = dicom.read_file(dicom_path)
                try:
                    orien_mat[index, ...] = dcm_header[0x0020,0x0032].value
                except:
                    LOGGER.warn('''dcm2nii -- {scan} -- {dicom} file could not be read properly (no (0020,0032) tag).'''.format(scan=scan_info['scan_id'],
                                                                                                                                dicom=dicom_path))
            # Get the axis with the biggest variance = axis used by the scanner
            var_mat = [np.var(orien_mat[:,0]), np.var(orien_mat[:,1]), np.var(orien_mat[:,2])]
            if var_mat[0] == np.max(var_mat):
                col = 0
            elif var_mat[1] == np.max(var_mat):
                col = 1
            else:
                col = 2

            # sort the matrix following the scanner axis
            sorted_orien_mat = orien_mat[orien_mat[:,col].argsort()]
            # compute the matrix of spacing between slices
            diff = sorted_orien_mat[:-1]-sorted_orien_mat[1:]
            spacing_mat = np.sqrt(np.multiply(diff, diff))

            #if this difference is bigger than 0.001, it means that one spacing
            #between two slices is bigger than the spacing between the closest
            #slices ---> some slices are missing
            max_array = [np.max(spacing_mat[:,0]), np.max(spacing_mat[:,1]), np.max(spacing_mat[:,2])]
            min_array = [np.min(spacing_mat[:,0]), np.min(spacing_mat[:,1]), np.min(spacing_mat[:,2])]
            max_spacing_mat = np.subtract(max_array, min_array)
            if max_spacing_mat[col] > 0.001:
                LOGGER.warn('''dcm2nii -- {scan} -- Slices might be missing or duplicated slices'''.format(scan=scan_info['scan_id']))
                self.log_warning_error('Slices might be missing or duplicated slices', scan_info)
                return False

        return True

    def dcm2nii(self, dicom_path):
        """ convert dicom to nifti using dcm2nii """
        LOGGER.debug('convert dcm to nii...')
        dcm2nii_cmd = '''{dcm2nii} -a n -e n -d n -g y -f n -n y -p n -v y -x n -r n {dicom}'''.format(dcm2nii=os.path.join(self.dcm2niipath, 'dcm2nii'), dicom=dicom_path)
        try:
            _ = sb.check_output(dcm2nii_cmd.split())
        except sb.CalledProcessError:
            LOGGER.debug("DCM --> NII conversion failed")
            return False

        return True

    def dcmdjpeg(self, dcm_dir):
        """ converting the dicom to jpeg dicoms """
        LOGGER.debug('run dcmdjpeg on the DICOMs.')
        for number, dicoms in enumerate(os.listdir(dcm_dir)):
            dcmdjpeg_cmd = '''{dcmdjpeg} {original_dcm} {new_dcm}'''.format(dcmdjpeg=os.path.join(self.dcmdjpegpath, 'dcmdjpeg'),
                                                                            original_dcm=os.path.join(dcm_dir, dicoms),
                                                                            new_dcm=os.path.join(dcm_dir, 'final_'+str(number)+'.dcm'))
            os.system(dcmdjpeg_cmd)
            os.remove(os.path.join(dcm_dir, dicoms))

    def upload_converted_images(self, dcm_dir, scan_obj, scan_info):
        """ upload the images after checking them """
        #Local variables
        nifti_list = []
        bval_fpath = ''
        bvec_fpath = ''

        #Get the bvec/bval files and NIFTI from the folder:
        for fpath in glob.glob(os.path.join(dcm_dir, '*')):
            if os.path.isfile(fpath):
                if fpath.lower().endswith('.bval'):
                    bval_fpath = fpath
                if fpath.lower().endswith('.bvec'):
                    bvec_fpath = fpath
                if fpath.lower().endswith('.nii.gz'):
                    nifti_list.append(fpath)
                if fpath.lower().endswith('.nii'):
                    os.system('gzip '+fpath)
                    nifti_list.append(fpath+'.gz')

        #Check NIFTI:
        good_to_upload = self.check_outputs(scan_info, nifti_list, bval_fpath, bvec_fpath)
        #Upload files:
        if good_to_upload:
            if os.path.isfile(bval_fpath) and os.path.isfile(bvec_fpath):
                #BVAL/BVEC
                XnatUtils.upload_file_to_obj(bval_fpath, scan_obj.resource('BVAL'), remove=True)
                XnatUtils.upload_file_to_obj(bvec_fpath, scan_obj.resource('BVEC'), remove=True)
                #keep the NII with the same name than the BVAL/BVEC
                nifti_list = filter(lambda x: x[:-7] == bval_fpath[:-5],nifti_list)
                XnatUtils.upload_files_to_obj(nifti_list, scan_obj.resource('NIFTI'), remove=True)
            else:
                #NII
                XnatUtils.upload_files_to_obj(nifti_list, scan_obj.resource('NIFTI'), remove=True)
            #more than one NIFTI uploaded
            if len(nifti_list) > 1:
                LOGGER.warn('''dcm2nii -- {scan} -- more than one NIFTI upload'''.format(scan=scan_info['scan_id']))
                self.log_warning_error('more than one NIFTI upload', scan_info)

    def check_outputs(self, scan_info, nifti_list, bval, bvec):
        """ Check that the outputs are right (opening nifti works)"""
        for nifti_fpath in nifti_list:
            try:
                nii = nib.load(nifti_fpath)
            except:
                LOGGER.warn('''dcm2nii -- {scan} -- {file} is not a proper NIFTI'''.format(scan=scan_info['scan_id'], file=os.path.basename(nifti_fpath)))
                self.log_warning_error('check the scan DICOM. No upload. a non-valid nifti was created', scan_info, error=True)
                return False
        if 'dti' in scan_info['type'].lower() or 'dif' in scan_info['type'].lower() or 'hardi' in scan_info['type'].lower():
            if not os.path.isfile(bval) or not os.path.isfile(bvec):
                LOGGER.warn('''dcm2nii -- {scan} -- no bval/bvec generated (DTI scan)'''.format(scan=scan_info['scan_id']))
                self.log_warning_error('check the scan DICOM. No upload. no bval/bvec were generated', scan_info, error=True)
                return False
        return True
