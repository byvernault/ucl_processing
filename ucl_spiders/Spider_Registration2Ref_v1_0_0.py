"""Spider to perform Diffusion Model Fitting.

Author:         Benjamin Yvernault
contact:        b.yvernault@ucl.ac.uk
Spider name:    Registration2Ref
Spider version: 1.0.0
Creation date:  2016-02-22 16:46:39.832273
Purpose:        Register chosen scans from a XNAT session to one of the scan
                chosen as a reference
"""

# Python packages import
import os
import sys
import glob
import numpy as np
import nibabel as nib
import subprocess as sb
from datetime import datetime
import matplotlib.pyplot as plt
from dax import spiders, XnatUtils, SessionSpider

__author__ = "Benjamin Yvernault"
__email__ = "b.yvernault@ucl.ac.uk"
__purpose__ = "Register chosen scans from a XNAT session to \
one of the scan chosen as a reference"
__spider_name__ = "Registration2Ref"
__version__ = "1.0.0"
__modifications__ = """2016-02-22 16:46:39.832273 - Original write
2016-05-10 18:04:01 - Update to new format respecting pep8
"""

REG_ALADIN_CMD = """{exe_path} \
--flo {input} \
--ref {reference} \
--aff {affine} \
--res {output} \
{args}
"""
RA_ARGS = """--floUpThr 1000\
--smooF 0 \
--maxit 5 \
--refUpThr 1000 \
--interp 0 \
--lp 3 \
--smooR 0 \
--ln 3 \
--pi 50 \
--refLowThr 0 \
--smooF 0 \
--maxit 5 \
--refUpThr 1000 \
--interp 0 \
--lp 3 \
--smooR 0 \
--ln 3 \
--pi 50 \
--refLowThr 0 \
--pv 50 \
--floLowThr 0
"""


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
        --scansType : scan types to register to the reference
        --scanRef : scan ID for the reference from XNAT
        --regAladin : path to reg_aladin's executable
        --openmp_core : number of core use by reg_aladin. Default: one

    :return: argument parser object created by parse_args()
    """
    ap = spiders.get_session_argparser("Registration2Ref", __purpose__)
    ap.add_argument("--scansID", dest="scans_id", required=True,
                    help="Scans ID from XNAT to register to reference.")
    ap.add_argument("--scanRef", dest="scan_ref", required=True,
                    help="Scan ID from XNAT of the reference \
scan for the registration.")
    ap.add_argument("--regAladin", dest="reg_aladin_exe", required=True,
                    help="path to reg_aladin's executable.")
    ap.add_argument("--openmp_core", dest="openmp_core", default=1,
                    help="Number of core used by reg_aladin.")
    return ap.parse_args()


class Spider_Registration2Ref(SessionSpider):
    """Register chosen scans to one of the scan chosen as a reference.

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
                 scan_ref, scans_id,
                 xnat_host=None, xnat_user=None, xnat_pass=None, suffix=""):
        """Entry point for Spider_Registration2Ref Class."""
        super(Spider_Registration2Ref, self).__init__(spider_path, jobdir,
                                                      xnat_project,
                                                      xnat_subject,
                                                      xnat_session,
                                                      xnat_host,
                                                      xnat_user,
                                                      xnat_pass,
                                                      suffix)
        self.xnat_scan_ref = scan_ref
        self.scans_id = scans_id
        self.inputs = list()
        self.reference = list()
        self.outputs = list()
        self.pdf_final = os.path.join(self.jobdir, 'Registration2Ref.pdf')
        # Print version for Niftyreg - GIFi
        pversion = sb.Popen([ARGS.reg_aladin_exe, '--version'],
                            stdout=sb.PIPE,
                            stderr=sb.PIPE)
        nve_version, _ = pversion.communicate()
        self.time_writer('reg_aladin (Niftyreg) version: %s' %
                         (nve_version.strip()))

    def pre_run(self):
        """Method to download data from XNAT."""
        resource = 'NIFTI'  # resource to download from the scan on XNAT
        input_folder = os.path.join(self.jobdir, 'inputs')
        if not os.path.exists(input_folder):
            os.makedirs(input_folder)

        for scan_id in self.scans_id:
            scan_path = os.path.join(input_folder, scan_id)
            if not os.path.exists(scan_path):
                os.makedirs(scan_path)
            self.inputs.extend(self.download(scan_id, resource, scan_path))

        # Find and select reference:
        self.reference.extend(self.download(self.xnat_scan_ref, resource,
                                            input_folder))

    def run(self, reg_aladin_exe):
        """Method running the process for the spider on the inputs data."""
        if reg_aladin_exe.endswith('reg_aladin'):
            exe_path = reg_aladin_exe
        elif os.path.isdir(reg_aladin_exe):
            exe_path = os.path.join(reg_aladin_exe, "reg_aladin")

        if not os.path.exists(exe_path):
            raise Exception("Executable '%s' not found" % (exe_path))
        else:
            output_folder = os.path.join(self.jobdir, 'outputs')
            if not os.path.exists(output_folder):
                os.makedirs(output_folder)

            for scan_id in self.scans_id:
                prefix_str = '%s_reg2_%s' % (scan_id, self.xnat_scan_ref)
                reg_folder = os.path.join(output_folder, prefix_str)
                if not os.path.exists(reg_folder):
                    os.makedirs(reg_folder)
                ipath = os.path.join(self.jobdir, 'inputs', scan_id, 'NIFTI',
                                     '*.nii.gz')
                input_fpath = glob.glob(ipath)[0]
                self.time_writer("Register %s to %s" %
                                 (input_fpath, self.reference[0]))
                output_name = "%s.nii" % (prefix_str)
                affine_name = "%s_affine_transformation.txt" % (prefix_str)
                output_fpath = os.path.join(reg_folder, output_name)
                affine_fpath = os.path.join(reg_folder, affine_name)
                cmd = REG_ALADIN_CMD.format(exe_path=exe_path,
                                            reference=self.reference[0],
                                            input=input_fpath,
                                            output=output_fpath,
                                            affine=affine_fpath,
                                            args=RA_ARGS)
                self.run_system_cmd(cmd)
                XnatUtils.gzip_nii(reg_folder)
                self.outputs.append(reg_folder)
            self.make_pdf()

    def finish(self):
        """Method to copy the results in dax.RESULTS_DIR."""
        results_dict = {'PDF': self.pdf_final}
        for out_folder in self.outputs:
            results_dict[os.path.basename(out_folder)] = out_folder
        self.upload_dict(results_dict)
        self.end()

    def plot_images(self, fig, data, number_images, subplot_index, label):
        """Method to plot the images on the PDF for the first page.

        :param fig: figure from matplotlib
        :param data: numpy array of the images (3D) to display
        :param number_images: number of images for the subplot (6)
        :param subplot_index: index of the line to display (1-6)
        :param label: y-label
        :return: None
        """
        ax = fig.add_subplot(number_images, 3, 3*subplot_index+1)
        ax.imshow(np.rot90(data[:, :, data.shape[2]/2]), cmap='gray')
        if subplot_index == 0:
            ax.set_title('Axial', fontsize=7)
        ax.set_ylabel(label, fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
        ax = fig.add_subplot(number_images, 3, 3*subplot_index+2)
        ax.imshow(np.rot90(data[:, data.shape[1]/2, :]), cmap='gray')
        if subplot_index == 0:
            ax.set_title('Coronal', fontsize=7)
        ax.set_axis_off()
        ax = fig.add_subplot(number_images, 3, 3*subplot_index+3)
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
        plt.figtext(0.5, 0.985, '-- Registration2Ref Pipeline PDF report --',
                    horizontalalignment='center', fontsize=10)
        footer = 'Date: %s -- page %s/%s -- PDF generated by TIG \
laboratory at UCL, London' % (str(date), str(pdf_page_number), str(pdf_pages))
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
        date = datetime.now()
        nb_images_per_page = 5

        # PDF path:
        pdf_pages_list = list()

        pdf_page_number = 0  # 5 images per pages
        fig = plt.figure(0, figsize=(7.5, 10))
        total_nb_pages = (len(self.outputs) + 1)/nb_images_per_page
        if total_nb_pages == 0:
            total_nb_pages = 1
        # Plot the reference
        # Open niftis with nibabel
        f_img = nib.load(self.reference[0])
        f_img_data = f_img.get_data()
        # Draw
        if len(f_img_data.shape) == 3:
            data = f_img_data
        elif len(f_img_data.shape) == 4:
            data = f_img_data[:, :, :, f_img_data.shape[3]/2]
        self.plot_images(fig, data, nb_images_per_page, 0, 'Reference')

        page_0_saved = False
        for index, out_folder in enumerate(self.outputs):
            # Get the output images
            nii_file = glob.glob(os.path.join(out_folder, '*.nii.gz'))
            if len(nii_file) != 1:
                err = '%s output image not found or more than one file found.'
                raise Exception(err % (nii_file))

            if (index+1)/5 != pdf_page_number:
                page_0_saved = True
                pdf_name = 'Registration2Ref_page%s.pdf' % (pdf_page_number)
                pdf_page = os.path.join(self.jobdir, pdf_name)
                # Save figure to PDF page
                self.save_pdf_page(fig, pdf_page, date, pdf_page_number,
                                   total_nb_pages)
                pdf_pages_list.append(pdf_page)
                pdf_page_number += 1
                fig = plt.figure(pdf_page_number, figsize=(7.5, 10))

            # Open niftis with nibabel
            f_img = nib.load(nii_file[0])
            f_img_data = f_img.get_data()
            # Draw
            if len(f_img_data.shape) == 3:
                data = f_img_data
            elif len(f_img_data.shape) == 4:
                data = f_img_data[:, :, :, f_img_data.shape[3]/2]
            self.plot_images(fig, data, nb_images_per_page, index+1,
                             os.path.basename(nii_file[0]).split('.')[0])

        # If only 5 or less images, save PDF zero
        if not page_0_saved:
            self.save_pdf_page(fig, self.pdf_final, date, 1, total_nb_pages)
        else:
            # Join the two pages for the PDF:
            cmd = 'gs -q -sPAPERSIZE=letter -dNOPAUSE -dBATCH \
-sDEVICE=pdfwrite -sOutputFile=%s' % (self.pdf_final)
            for page in pdf_pages_list:
                cmd = '%s %s' % (cmd, page)
            self.time_writer('INFO:saving final PDF: %s ' % cmd)
            os.system(cmd)

if __name__ == '__main__':
    ARGS = parse_args()
    # Set the environment variable:
    os.environ['OMP_NUM_THREADS'] = ARGS.openmp_core

    # generate spider object:
    spider_obj = Spider_Registration2Ref(spider_path=sys.argv[0],
                                         jobdir=ARGS.temp_dir,
                                         xnat_project=ARGS.proj_label,
                                         xnat_subject=ARGS.subj_label,
                                         xnat_session=ARGS.sess_label,
                                         scan_ref=ARGS.scan_ref,
                                         scans_id=ARGS.scans_id.split(','),
                                         xnat_host=ARGS.host,
                                         xnat_user=ARGS.user,
                                         xnat_pass=None,
                                         suffix=ARGS.suffix)
    # print some information before starting
    spider_obj.print_init(ARGS, "Benjamin Yvernault", "b.yvernault@ucl.ac.uk")

    # Pre-run method to download data from XNAT
    spider_obj.pre_run()

    # Run method
    spider_obj.run(ARGS.reg_aladin_exe)

    # Finish method to copy results
    spider_obj.finish()
