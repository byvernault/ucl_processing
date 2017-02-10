"""Spider to execute Vessel Extraction (UCL).

Author:         Benjamin Yvernault
contact:        b.yvernault@ucl.ac.uk
Spider name:    Vessel_Extraction
Spider version: 1.0.0
Creation date:  2016-02-22 14:19:24.698923
Purpose:        Extract the vessel from the T1/MPRAGE scan
"""

import os
import sys
import numpy as np
import nibabel as nib
from datetime import datetime
import matplotlib.pyplot as plt
from dax import XnatUtils, spiders, ScanSpider

__author__ = "Benjamin Yvernault"
__email__ = "b.yvernault@ucl.ac.uk"
__purpose__ = "Extract the vessel from the T1/MPRAGE scan"
__spider_name__ = "Vessel_Extraction"
__version__ = "1.0.0"
__modifications__ = """2016-02-22 14:19:24.698923 - Original write
2016-05-10 18:04:01 - Update to new format respecting pep8
"""

DEFAULT_PIXEL_SIZE = '0.775438'
DEFAULT_EXE_PATH = 'niftkVesselExtractor'
VESSEL_CMD = """{exe_path} -i {input} -o {output} --mod 0 --aone 0.5 --atwo 2 \
--min {pixel_size} --max 3.09 --intfil"""


def parse_args():
    """Argument parser for the spider input variables.

    by default (set in get_session_argparser):
        -p       : proj_label
        -s       : subj_label
        -e       : sess_label
        -c       : scan_label
        -d       : temp_dir
        --suffix : suffix (suffix for assessor proctype on XNAT)
        --host   : host for XNAT (default: XNAT_HOST env variable)
        --user   : user on XNAT (default: XNAT_USER env variable)
    your arguments:
        --vesselExtPath : path to niftkVesselExtractor's executable

    :return: argument parser object created by parse_args()
    """
    ap = spiders.get_scan_argparser("Vessel_Extraction", __purpose__)
    ap.add_argument("--vesselExtPath", dest="vessel_extractor",  required=True,
                    help="path to the executable niftkVesselExtractor.")
    ap.add_argument("--drc", dest="drc_ana",  required=False,
                    action='store_true',
                    help="Set the NIFTK_DRC_ANALYZE to ON.")
    return ap.parse_args()


class Spider_Vessel_Extraction(ScanSpider):
    """Scan Spider class to do: Extract the vessel from the T1/MPRAGE scan.

    :param spider_path: spider file path
    :param jobdir: directory for temporary files
    :param xnat_project: project ID on XNAT
    :param xnat_subject: subject label on XNAT
    :param xnat_session: experiment label on XNAT
    :param xnat_scan: scan label on XNAT
    :param xnat_host: host for XNAT if not set in environment variables
    :param xnat_user: user for XNAT if not set in environment variables
    :param xnat_pass: password for XNAT if not set in environment variables
    :param exe_path: niftkVesselExtractor path
    :param suffix: suffix to the assessor creation
    """

    def __init__(self, spider_path, jobdir, xnat_project, xnat_subject,
                 xnat_session, xnat_scan, xnat_host=None, xnat_user=None,
                 xnat_pass=None, exe_path=DEFAULT_EXE_PATH, suffix=""):
        """Entry point for Spider_Vessel_Extraction Class."""
        super(Spider_Vessel_Extraction, self).__init__(spider_path,
                                                       jobdir,
                                                       xnat_project,
                                                       xnat_subject,
                                                       xnat_session,
                                                       xnat_scan,
                                                       xnat_host,
                                                       xnat_user,
                                                       xnat_pass,
                                                       suffix)
        self.inputs = list()
        self.output = ''
        self.pdfpath = ''
        self.pixel_size = DEFAULT_PIXEL_SIZE
        self.exe_path = self.check_executable(exe_path, 'niftkVesselExtractor')

    def get_voxel_size(self):
        """Method to get the voxel size from XNAT if define using default value if not.

        :return: voxel size [x, y, z]
        """
        xnat = XnatUtils.get_interface(host=self.host,
                                       user=self.user,
                                       pwd=self.pwd)
        scan_obj = self.select_obj(xnat, self.xnat_scan, None)
        vsize = scan_obj.attrs.mget(['xnat:mrScanData/parameters/voxelRes/x',
                                     'xnat:mrScanData/parameters/voxelRes/y',
                                     'xnat:mrScanData/parameters/voxelRes/z'])
        xnat.disconnect()
        return vsize

    def pre_run(self):
        """Method to download data from XNAT.

        :return: None
        """
        resource = 'NIFTI'  # resource to download from the scan on XNAT
        folder = os.path.join(self.jobdir, 'inputs')
        if not os.path.exists(os.path.join(self.jobdir, 'inputs')):
            os.makedirs(folder)
        self.inputs.extend(self.download(self.xnat_scan, resource, folder))

        # Read from XNAT the pixel size from DICOM, if not default value:
        vsize = self.get_voxel_size()
        if not vsize:
            msg = "Using default value for pixel size in niftkVesselExtractor.\
Value not set on XNAT."
            self.time_writer(msg)
        else:
            self.pixel_size = min(vsize)

    def run(self):
        """Method running the process for the spider on the inputs data.

        :return: None
        """
        if not os.path.exists(os.path.join(self.jobdir, 'outputs')):
            os.makedirs(os.path.join(self.jobdir, 'outputs'))
        output_name = ("%s_vessel_extracted.nii.gz" %
                       (os.path.basename(self.inputs[0]).split('.')[0]))
        self.output = os.path.join(self.jobdir, 'outputs', output_name)
        cmd = VESSEL_CMD.format(exe_path=self.exe_path,
                                input=self.inputs[0],
                                output=self.output,
                                pixel_size=self.pixel_size)
        self.run_system_cmd(cmd)
        self.make_pdf()

    def finish(self):
        """Method to copy the results in dax.RESULTS_DIR.

        :return: None
        """
        results_dict = {
                        'PDF': self.pdfpath,
                        'OUTPUT': self.output
                        }
        self.upload_dict(results_dict)
        self.end()

    def make_pdf(self):
        """Method to make pdf."""
        f_img = nib.load(self.output)
        data = f_img.get_data()
        nb_slices = data.shape[2]
        slices = {'0': [(nb_slices-1)/12, (nb_slices-1)/6,
                        (nb_slices-1)/4, (nb_slices-1)/3],
                  '1': [5*(nb_slices-1)/12, (nb_slices-1)/2,
                        7*(nb_slices-1)/12, 2*(nb_slices-1)/3],
                  '2': [9*(nb_slices-1)/12, 5*(nb_slices-1)/6,
                        11*(nb_slices-1)/12, (nb_slices-1)]}
        labels = {'0': '',
                  '1': '',
                  '2': ''}
        vmins = {'0': 0,
                 '1': 0,
                 '2': 0}
        vmaxs = {'0': 100,
                 '1': 100,
                 '2': 100}
        images = [self.output, self.output, self.output]
        self.pdfpath = os.path.join(self.jobdir, 'vessel_report.pdf')
        self.plot_images_page(self.pdfpath, 1, images, 'Vessel Extraction',
                              labels, slices=slices, cmap='hot',
                              vmins=vmins, vmaxs=vmaxs)

    def plot_images_page(self, pdf_path, page_index, nii_images, title,
                         image_labels, slices=None, cmap='gray',
                         vmins=None, vmaxs=None):
        """Plot list of images (3D-4D) on a figure (PDF page).

        plot_images_figure will create one pdf page with only images.
        Each image corresponds to one line with by default axial/sag/cor view
        of the mid slice. If you use slices, it will show different slices of
        the axial plan view. You can specify the cmap and the vmins/vmaxs if
        needed by using a dictionary with the index of each line (0, 1, ...).

        :param pdf_path: path to the pdf to save this figure to
        :param page_index: page index for PDF
        :param nii_images: python list of nifty images
        :param title: Title for the report page
        :param image_labels: list of titles for each images
            one per image in nii_images
        :param slices: dictionary of list of slices to display
            if None, display axial, coronal, sagital
        :param cmap: cmap to use to display images or dict
            of cmaps for each images with the indices as key
        :param vmins: define vmin for display (dict)
        :param vmaxs: define vmax for display (dict)
        :return: pdf path created

        E.g for two images:
        images = [imag1, image2]
        slices = {'0':[50, 80, 100, 130],
                  '1':[150, 180, 200, 220]}
        labels = {'0': 'Label 1',
                  '1': 'Label 2'}
        cmaps = {'0':'hot',
                 '1': 'gray'}
        vmins = {'0':10,
                 '1':20}
        vmaxs = {'0':100,
                 '1':150}
        """
        self.time_writer('INFO: generating pdf page %d with images.'
                         % page_index)
        fig = plt.figure(page_index, figsize=(7.5, 10))
        # Titles:
        if not isinstance(cmap, dict):
            default_cmap = cmap
            cmap = {}
        else:
            default_cmap = 'gray'
        if not isinstance(vmins, dict):
            self.time_writer("Warning: vmins wasn't a dictionary. \
Using default.")
            vmins = {}
        if not isinstance(vmins, dict):
            self.time_writer("Warning: vmaxs wasnt' a dictionary. \
Using default.")
            vmaxs = {}
        if isinstance(nii_images, str):
            nii_images = [nii_images]
        number_im = len(nii_images)
        for index, image in enumerate(nii_images):
            # Open niftis with nibabel
            f_img = nib.load(image)
            data = f_img.get_data()
            if len(data.shape) == 4:
                data = data[:, :, :, data.shape[3]/2]
            default_slices = [data.shape[2]/4, data.shape[2]/2,
                              3*data.shape[2]/4]
            default_label = 'Line %s' % index

            if slices:
                if not isinstance(slices, dict):
                    self.time_writer("Warning: slices wasn't a dictionary. \
Using default.")
                    slices = {}
                self.time_writer('INFO: showing different slices.')
                li_slices = slices.get(str(index), default_slices)
                slices_number = len(li_slices)
                for slice_ind, slice_value in enumerate(li_slices):
                    ind = slices_number*index+slice_ind+1
                    ax = fig.add_subplot(number_im, slices_number, ind)
                    data_z_rot = np.rot90(data[:, :, slice_value])
                    ax.imshow(data_z_rot,
                              cmap=cmap.get(str(index), default_cmap),
                              vmin=vmins.get(str(index), None),
                              vmax=vmaxs.get(str(index), None))
                    ax.set_title('Slice %d' % slice_value, fontsize=7)
                    ax.set_xticks([])
                    ax.set_yticks([])
                    if slice_ind == 0:
                        ax.set_ylabel(image_labels.get(str(index),
                                      default_label), fontsize=9)
            else:
                self.time_writer('INFO: display different plan view \
(ax/sag/cor) of the mid slice.')
                ax = fig.add_subplot(number_im, 3, 3*index+1)
                data_z_rot = np.rot90(data[:, :, data.shape[2]/2])
                ax.imshow(data_z_rot, cmap=cmap.get(str(index), default_cmap),
                          vmin=vmins.get(str(index), None),
                          vmax=vmaxs.get(str(index), None))
                ax.set_title('Axial', fontsize=7)
                ax.set_ylabel(image_labels.get(str(index), default_label),
                              fontsize=9)
                ax.set_xticks([])
                ax.set_yticks([])
                ax = fig.add_subplot(number_im, 3, 3*index+2)
                data_y_rot = np.rot90(data[:, data.shape[1]/2, :])
                ax.imshow(data_y_rot, cmap=cmap.get(str(index), default_cmap),
                          vmin=vmins.get(str(index), None),
                          vmax=vmaxs.get(str(index), None))
                ax.set_title('Coronal', fontsize=7)
                ax.set_axis_off()
                ax = fig.add_subplot(number_im, 3, 3*index+3)
                data_x_rot = np.rot90(data[data.shape[0]/2, :, :])
                ax.imshow(data_x_rot, cmap=cmap.get(str(index), default_cmap),
                          vmin=vmins.get(str(index), None),
                          vmax=vmaxs.get(str(index), None))
                ax.set_title('Sagittal', fontsize=7)
                ax.set_axis_off()

        fig.tight_layout()
        date = datetime.now()
        # Titles page
        plt.figtext(0.5, 0.985, '-- %s PDF report --' % title,
                    horizontalalignment='center', fontsize=12)
        plt.figtext(0.5, 0.02, 'Date: %s -- page %d' % (str(date), page_index),
                    horizontalalignment='center', fontsize=8)
        plt.show()
        fig.savefig(pdf_path, transparent=True, orientation='portrait',
                    dpi=100)
        plt.close(fig)
        return pdf_path


if __name__ == '__main__':
    ARGS = parse_args()
    # Set environment variable for fixing flips
    if ARGS.drc_ana:
        os.environ["NIFTK_DRC_ANALYZE"] = 'ON'

    # generate spider object:
    spider_obj = Spider_Vessel_Extraction(spider_path=sys.argv[0],
                                          jobdir=ARGS.temp_dir,
                                          xnat_project=ARGS.proj_label,
                                          xnat_subject=ARGS.subj_label,
                                          xnat_session=ARGS.sess_label,
                                          xnat_scan=ARGS.scan_label,
                                          xnat_host=ARGS.host,
                                          xnat_user=ARGS.user,
                                          xnat_pass=None,
                                          exe_path=ARGS.vessel_extractor,
                                          suffix=ARGS.suffix)
    # print some information before starting
    spider_obj.print_init(ARGS, "Benjamin Yvernault", "b.yvernault@ucl.ac.uk")

    # Pre-run method to download data from XNAT
    spider_obj.pre_run()

    # Check the value:
    spider_obj.time_writer("NIFTK_DRC_ANALYZE set to %s"
                           % os.environ.get("NIFTK_DRC_ANALYZE", 'OFF'))

    # Run method
    spider_obj.run()

    # Finish method to copy results
    if not ARGS.skipfinish:
        spider_obj.finish()
