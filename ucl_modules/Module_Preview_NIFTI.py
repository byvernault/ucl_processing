""" Module to generate the Preview for a scan from a resource ( hosting NIFTI file) """
from dax import XnatUtils, ScanModule
import os
import logging
import scipy as sp
import numpy as np
import nibabel as nib

LOGGER = logging.getLogger('dax')

DEFAULT_TPM_PATH = '/tmp/Preview_NIFTI_tmp/'
DEFAULT_MODULE_NAME = 'Preview_NIFTI'
DEFAULT_TEXT_REPORT = 'ERROR/WARNING for preview NIFTI :\n'
DEFAULT_RESOURCE = 'NIFTI'

class Module_Preview_NIFTI(ScanModule):
    """ Module to generate preview for XNAT """
    def __init__(self, mod_name=DEFAULT_MODULE_NAME, directory=DEFAULT_TPM_PATH, email=None, text_report=DEFAULT_TEXT_REPORT, resourcename=DEFAULT_RESOURCE):
        """ init function overridden from base-class"""
        super(Module_Preview_NIFTI, self).__init__(mod_name, directory, email, text_report=text_report)
        if isinstance(resourcename, str):
            self.resourcename = resourcename
        else:
            self.resourcename = DEFAULT_RESOURCE

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
            LOGGER.warn('Preview NIFTI -- afterrun -- '+self.directory+' not empty. Could not delete it.')

    def needs_run(self, cscan, xnat):
        """ needs_run function overridden from base-class
                cscan = CacheScan object from XnatUtils
            return True or False
        """
        scan_info = cscan.info()
        if XnatUtils.has_resource(cscan, 'SNAPSHOTS'):
            LOGGER.debug('Has snapshot')
            return False

        if not XnatUtils.has_resource(cscan, self.resourcename):
            LOGGER.warn('Preview NIFTI -- '+scan_info['scan_id']+' -- no '+self.resourcename+' resource')
            return False

        return True

    def run(self, scan_info, scan_obj):
        """ main function to run on scan """
        if not len(scan_obj.resource(self.resourcename).files().get()) > 0:
            LOGGER.warn('Preview NIFTI -- '+scan_info['scan_id']+' -- No file(s) for '+self.resourcename)
        else:
            LOGGER.debug("downloading "+self.resourcename+"...")

            fpath = XnatUtils.download_biggest_file_from_obj(self.directory, scan_obj.resource(self.resourcename))

            if not os.path.exists(fpath):
                LOGGER.warn('Preview NIFTI -- '+scan_info['scan_id']+' -- '+self.resourcename+' file size is zero.')
            else:
                sm_gif = str(os.path.join(self.directory, os.path.splitext(os.path.basename(fpath))[0] + "_sm.gif"))
                lg_gif = str(os.path.join(self.directory, os.path.splitext(os.path.basename(fpath))[0] + "_lg.gif"))

                # Generate preview from the nifti file
                self.generate_preview(fpath, sm_gif, lg_gif)

                if os.path.isfile(sm_gif):
                    LOGGER.debug('GIF Made / upload ...')
                    scan_obj.resource('SNAPSHOTS').file('snap_t.gif').put(sm_gif, 'GIF', 'THUMBNAIL')
                    scan_obj.resource('SNAPSHOTS').file('snap.gif').put(lg_gif, 'GIF', 'ORIGINAL')
                else:
                    LOGGER.error('Preview NIFTI -- '+scan_info['scan_id']+' -- GIF FAILED')
                    self.log_warning_error('GIF failed for NIFTI', scan_info, error=True)

                self.clean_directory()

    def generate_preview(self, nifti, smgif, lggif):
        """
        generate snapshots of the nifti given as parameters
        :param nifti: path to the nifti files
        :return: path to the original and preview SNAPSHOTs of the NIFTI
        """
        print "   --> Generating snapshots "+nifti
        global out
        nii = self.load_nifti(nifti)

        if nii:
            # fix direction from the nifti to have proper snapshots
            data = self.get_data(nii)
            size = data.shape
            dtype = str(data.dtype)
            if size[2] > 100:
                # Keep 100 slices from the nifti
                data = self.resize_100_slices(data, size)
            # Expanding the matrix
            data = np.expand_dims(data, axis=3)
            # transpose third and fourth dimension
            final_data = np.transpose(data, (0, 1, 3, 2))
            if "float" in dtype.lower():
                # rescale our final data
                final_data = self.rescale(final_data)
            # Generate our big size1 * size2 numpyArray for snapshots
            size = final_data.shape
            n2 = np.ceil(np.sqrt(size[3]))
            out = np.zeros([n2*size[0], n2*size[1]])
            row = 0
            col = 0
            for i in range(size[3]):
                if col >= n2:
                    row += 1
                    col = 0

                out[(size[0]*row):(size[0]+size[0]*row), (size[1]*col):(size[1]+size[1]*col)] = final_data[:, :, 0, i]
                col += 1
        # Save the images as gif for small and large snapshots:
        sp.misc.imsave(smgif, sp.misc.imresize(np.multiply(out, 255), (512, 512)))
        sp.misc.imsave(lggif, np.multiply(out, 255))

    @staticmethod
    def get_data(nii_obj):
        """
        Extract the data from the nifti.
        Rotate the images depending on the orientation
        :param nifti: nibabel nifti object
        :return: return data from the nifti
        """
        data = nii_obj.get_data()
        # rotating if bitpix is 16 (hack)
        if nii_obj.header['bitpix'] == 16:
            print "   --> Rotating 90 degres slices"
            data = np.rot90(data)
        # Reshape if 4D images
        size = data.shape
        if len(size) == 4:
            # Reshape the numpay array in a 3D array
            # (swap axes to keep matlab standard)
            sl = reduce(lambda x, y: x*y, size[2:])
            data = np.swapaxes(data, 2, 3).reshape((size[0], size[0], sl))
        return data

    @staticmethod
    def load_nifti(nifti):
        """
        Load nifti as a nibabel object
        :param nifti: path to the nifti files
        :return: return NIFTI object or NULL
        """
        if not os.path.exists(nifti):
            print "load_nifti - error - couldn't find the file "+nifti
            return None
        try:
            nii = nib.load(nifti)
            return nii
        except:
            print "load_nifti - error - couldn't open NIFTI "+nifti
            return None

    @staticmethod
    def resize_100_slices(data, size):
        """
        Resize array to keep only 100 slices.
        :param data: numpy array of the data
        :param size: size of the data [x, y, z]
        :return: new array
        """
        print "   --> Limiting to 100 Slices"
        slices_list = np.round(np.linspace(0, size[2]-1, num=100))
        new_size = [size[0], size[1], len(slices_list)]
        new_array = np.zeros(new_size)
        #looping through array
        for index, i in enumerate(slices_list):
            new_array[..., index] = data[..., i]

        return new_array

    @staticmethod
    def rescale(data):
        """
        rescale python Array with user-defined min and max
        :param data: numpy array
        :return: array rescale
        """
        # Min and Max
        a_max = np.max(data)
        a_min = np.min(data)
        # Rescale the data between 0.0 - 1.0
        rescale_array = (data-a_min)/(a_max/2-a_min)
        rescale_array[rescale_array > 1.0] = 1.0

        return rescale_array
