"""Spider to perform gif parcellation.

Author:         Benjamin Yvernault
contact:        b.yvernault@ucl.ac.uk
Spider name:    GIF_Parcellation
Spider version: 1.0.0
Creation date:  2016-02-05 13:52:08.737570
Purpose:        Parcellation of the brain using GIF: Geodesic Information Flow.
"""

# Python packages import
import os
import sys
import csv
import glob
import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt
import subprocess as sb
from datetime import datetime
from dax import spiders, ScanSpider

__author__ = "Benjamin Yvernault"
__email__ = "b.yvernault@ucl.ac.uk"
__purpose__ = "Parcellation of the brain using GIF: Geodesic Information Flow."
__spider_name__ = "GIF_Parcellation"
__version__ = "1.0.0"
__modifications__ = """2016-03-15 14:56 - Adding working_dir options
2016-05-10 18:04:01 - Update to new format respecting pep8
"""


YLABELS = ['Bias Corrected', 'Brain', 'labels', 'Segmentation', 'tiv', 'prior']
CMAPS = ['gray', 'gray', None, 'gray', 'gray', None]
GIF_CMD = """python {exe_path} \
-i {input} \
-o {output} \
-d {db_xml} \
--no_qsub \
--openmp_core {number_core} \
--n_procs 1 \
--working_dir '{working_dir}'"""


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
        --dbt    : gif-based database xml file describing the inputs
        --gif    : Path to the Gif python script: perform_gif_propagation.py
        --openmp_core : number of core use by reg_aladin. Default: one
        --working_dir : working directory for temp files
    :return: argument parser object created by parse_args()
    """
    ap = spiders.get_scan_argparser("GIF_Parcellation", __purpose__)
    ap.add_argument("--dbt", dest="dbtemplate", required=True,
                    help="gif-based database xml file describing the inputs.")
    ap.add_argument("--gif", dest="gif_script", required=True,
                    help="Path to perform_gif_propagation.py.")
    ap.add_argument("--openmp_core", dest="openmp_core", default=1,
                    help="Number of core used by reg_aladin.")
    ap.add_argument("--working_dir", dest="working_dir", default="",
                    help="working directory for temp files. Default: output")
    return ap.parse_args()


class Spider_GIF_Parcellation(ScanSpider):
    """Spider for Parcellation of the brain using Geodesic Information Flow.

    :param spider_path: spider file path
    :param jobdir: directory for temporary files
    :param xnat_project: project ID on XNAT
    :param xnat_subject: subject label on XNAT
    :param xnat_session: experiment label on XNAT
    :param xnat_scan: scan label on XNAT
    :param xnat_host: host for XNAT if not set in environment variables
    :param xnat_user: user for XNAT if not set in environment variables
    :param xnat_pass: password for XNAT if not set in environment variables
    :param number_core: number of core used by reg_aladin
    :param suffix: suffix to the assessor creation
    """

    def __init__(self, spider_path, jobdir, xnat_project,
                 xnat_subject, xnat_session, xnat_scan,
                 xnat_host=None, xnat_user=None, xnat_pass=None,
                 number_core=1, suffix=""):
        """Entry point for Spider_GIF_Parcellation Class."""
        super(Spider_GIF_Parcellation, self).__init__(spider_path, jobdir,
                                                      xnat_project,
                                                      xnat_subject,
                                                      xnat_session,
                                                      xnat_scan,
                                                      xnat_host,
                                                      xnat_user,
                                                      xnat_pass,
                                                      suffix, subdir=False)
        self.inputs = list()
        self.number_core = number_core
        # Print version for Niftyreg - GIFi
        proc_nifty = sb.Popen(['reg_aladin', '--version'],
                              stdout=sb.PIPE,
                              stderr=sb.PIPE)
        nifty_reg_version, _ = proc_nifty.communicate()
        self.time_writer('Niftyreg version: %s' % (nifty_reg_version.strip()))
        self.time_writer('GIF version: 3a76a255ab')

    def pre_run(self):
        """Method to download data from XNAT.

        :return: None
        """
        resource = 'NIFTI'  # resource to download from the scan on XNAT
        folder = os.path.join(self.jobdir, 'inputs')
        if not os.path.exists(os.path.join(self.jobdir, 'inputs')):
            os.makedirs(folder)
        self.inputs.extend(self.download(self.xnat_scan, resource, folder))

    def run(self, gif_script, dbtemplate, working_dir):
        """Method running the process for the spider on the inputs data.

        :param gif_script: Path to the Gif python script
        :param dbtemplate: gif-based database xml file describing the inputs
        :return: None
        """
        if len(working_dir) > 0 and not os.path.exists(working_dir):
            os.makedirs(working_dir)

        if gif_script.endswith('perform_gif_propagation.py'):
            exe_path = gif_script
        else:
            exe_path = os.path.join(gif_script, "perform_gif_propagation.py")

        if not os.path.exists(exe_path):
            raise Exception("Python script: %s not found" % (exe_path))
        else:
            if not os.path.exists(os.path.join(self.jobdir, 'outputs')):
                os.makedirs(os.path.join(self.jobdir, 'outputs'))
            cmd = GIF_CMD.format(exe_path=exe_path,
                                 input=self.inputs[0],
                                 output=os.path.join(self.jobdir, 'outputs'),
                                 db_xml=dbtemplate,
                                 number_core=self.number_core,
                                 working_dir=working_dir)
            self.run_system_cmd(cmd)
            self.make_pdf()

    def finish(self):
        """Method to copy the results in dax.RESULTS_DIR."""
        out_dir = os.path.join(self.jobdir, 'outputs')
        # Images outputs:
        bias_corrected = glob.glob(os.path.join(out_dir,
                                                '*bias_corrected.nii.gz'))
        brain = glob.glob(os.path.join(out_dir, '*brain.nii.gz'))
        labels = glob.glob(os.path.join(out_dir, '*labels.nii.gz'))
        prior = glob.glob(os.path.join(out_dir, '*prior.nii.gz'))
        seg = glob.glob(os.path.join(out_dir, '*seg.nii.gz'))
        tiv = glob.glob(os.path.join(out_dir, '*tiv.nii.gz'))
        # Volumes:
        volumes = glob.glob(os.path.join(out_dir, '*volumes.csv'))
        results_dict = {'PDF': os.path.join(self.jobdir,
                                            'GIF_parcellation.pdf'),
                        'BIAS_COR': bias_corrected,
                        'BRAIN': brain,
                        'LABELS': labels,
                        'PRIOR': prior,
                        'SEG': seg,
                        'TIV': tiv,
                        'VOLUME': volumes}
        self.upload_dict(results_dict)
        self.end()

    def plot_images(self, fig, data, number_images, subplot_index):
        """Method to plot the images on the PDF for the first page.

        :param fig: figure from matplotlib
        :param data: numpy array of the images (3D) to display
        :param number_images: number of images for the subplot (6)
        :param subplot_index: index of the line to display (1-6)
        :return: None
        """
        ax = fig.add_subplot(number_images, 3, 3*subplot_index+1)
        data_z_rot = np.rot90(data[:, :, data.shape[2]/2])
        ax.imshow(data_z_rot, cmap=CMAPS[subplot_index])
        if subplot_index == 0:
            ax.set_title('Axial', fontsize=7)
        ax.set_ylabel(YLABELS[subplot_index], fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
        ax = fig.add_subplot(number_images, 3, 3*subplot_index+2)
        data_y_rot = np.rot90(data[:, data.shape[1]/2, :])
        ax.imshow(data_y_rot, cmap=CMAPS[subplot_index])
        if subplot_index == 0:
            ax.set_title('Coronal', fontsize=7)
        ax.set_axis_off()
        ax = fig.add_subplot(number_images, 3, 3*subplot_index+3)
        data_x_rot = np.rot90(data[data.shape[0]/2, :, :])
        ax.imshow(data_x_rot, cmap=CMAPS[subplot_index])
        if subplot_index == 0:
            ax.set_title('Sagittal', fontsize=7)
        ax.set_axis_off()

    def make_pdf(self):
        """Method to make the PDF for the spider.

        :return: None
        """
        # Define output files
        # Variables:
        date = datetime.now()
        print self.jobdir
        # PDF path:
        pdf_page1 = os.path.join(self.jobdir, 'GIF_parcellation_page1.pdf')
        pdf_page2 = os.path.join(self.jobdir, 'GIF_parcellation_page2.pdf')
        pdf_final = os.path.join(self.jobdir, 'GIF_parcellation.pdf')

        # Images outputs:
        out_dir = os.path.join(self.jobdir, 'outputs')
        bias_corrected = glob.glob(os.path.join(out_dir,
                                                '*bias_corrected.nii.gz'))
        brain = glob.glob(os.path.join(out_dir, '*brain.nii.gz'))
        labels = glob.glob(os.path.join(out_dir, '*labels.nii.gz'))
        prior = glob.glob(os.path.join(out_dir, '*prior.nii.gz'))
        seg = glob.glob(os.path.join(out_dir, '*seg.nii.gz'))
        tiv = glob.glob(os.path.join(out_dir, '*tiv.nii.gz'))
        list_images = [bias_corrected, brain, labels, seg, tiv, prior]
        # Volumes:
        volumes = glob.glob(os.path.join(out_dir, '*volumes.csv'))

        # Page 1:
        fig = plt.figure(1, figsize=(7.5, 10))
        number_images = len(list_images)
        for index, image_file in enumerate(list_images):
            if len(image_file) != 1:
                err = '%s output image not found or more than one file found.'
                raise Exception(err % (image_file))
            # Open niftis with nibabel
            f_img = nib.load(image_file[0])
            f_img_data = f_img.get_data()
            # Draw
            if len(f_img_data.shape) == 3:
                data = f_img_data
            elif len(f_img_data.shape) == 4:
                data = f_img_data[:, :, :, f_img_data.shape[3]/2]
            self.plot_images(fig, data, number_images, index)

        # Set footer and title
        fig.tight_layout()
        plt.figtext(0.5, 0.985, '-- GIF_Parcellation Pipeline PDF report --',
                    horizontalalignment='center', fontsize=10)
        bottom_page = 'Date: %s -- page 1/2 -- PDF generated by \
TIG laboratory at UCL, London' % str(date)
        plt.figtext(0.5, 0.02, bottom_page, horizontalalignment='center',
                    fontsize=8)

        # Writing Figure
        plt.show()
        fig.savefig(pdf_page1,
                    transparent=True,
                    orientation='portrait',
                    dpi=100)
        plt.close(fig)

        # Page 2
        print volumes
        if len(volumes) != 1:
            err = '%s output csv file with information on volumes not found \
or more than one file found.'
            raise Exception(err % (volumes))
        with open(volumes[0], 'rb') as csvfileread:
            csvreader = csv.reader(csvfileread, delimiter=',')
            list_labels_name = csvreader.next()
            list_labels_volume = csvreader.next()

        cell_text = []
        for ind in range(len(list_labels_volume)):
            txt = smaller_str(list_labels_name[ind].strip().replace('"', ''),
                              size=30)
            cell_text.append([txt, "%.2f" % float(list_labels_volume[ind])])

        # Make the table
        fig = plt.figure(2, figsize=(7.5, 10))
        nb_vol = 47
        for i in range(3):
            ax = fig.add_subplot(1, 3, i+1)
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)
            ax.axis('off')
            the_table = ax.table(cellText=cell_text[nb_vol*i:nb_vol*(i+1)],
                                 colColours=[(0.8, 0.4, 0.4), (1.0, 1.0, 0.4)],
                                 colLabels=['Label name', 'Volume'],
                                 colWidths=[0.8, 0.32],
                                 loc='center',
                                 rowLoc='left',
                                 colLoc='left',
                                 cellLoc='left')

            the_table.auto_set_font_size(False)
            the_table.set_fontsize(6)

        # Set footer and title
        plt.figtext(0.5, 0.95,
                    '-- Volumes computed by GIF_Parcellation Pipeline --',
                    horizontalalignment='center', fontsize=10)
        bottom_page = 'Date: %s -- page 2/2 -- PDF generated by \
TIG laboratory at UCL, London' % str(date)
        plt.figtext(0.5, 0.02, bottom_page, horizontalalignment='center',
                    fontsize=8)

        # Writing Figure
        plt.show()
        fig.savefig(pdf_page2,
                    transparent=True,
                    orientation='portrait',
                    dpi=300)
        plt.close(fig)

        # Join the two pages for the PDF:
        cmd = 'gs -q -sPAPERSIZE=letter -dNOPAUSE -dBATCH \
-sDEVICE=pdfwrite -sOutputFile=%s %s %s' % (pdf_final, pdf_page1, pdf_page2)
        self.time_writer('INFO:saving final PDF: %s ' % cmd)
        os.system(cmd)


def smaller_str(str_option, size=10, end=False):
    """Method to shorten a string into a smaller size.

    :param str_option: string to shorten
    :param size: size of the string to keep (default: 10 characters)
    :param end: keep the end of the string visible (default beginning)
    :return: shortened string
    """
    if len(str_option) > size:
        if end:
            return '...%s' % (str_option[-size:])
        else:
            return '%s...' % (str_option[:size])
    else:
        return str_option


if __name__ == '__main__':
    # arguments
    args = parse_args()

    # generate spider object:
    spider_obj = Spider_GIF_Parcellation(spider_path=sys.argv[0],
                                         jobdir=args.temp_dir,
                                         xnat_project=args.proj_label,
                                         xnat_subject=args.subj_label,
                                         xnat_session=args.sess_label,
                                         xnat_scan=args.scan_label,
                                         xnat_host=args.host,
                                         xnat_user=args.user,
                                         xnat_pass=None,
                                         number_core=args.openmp_core,
                                         suffix=args.suffix)

    # print some information before starting
    spider_obj.print_init(args, "Benjamin Yvernault", "b.yvernault@ucl.ac.uk")

    # Pre-run method to download data from XNAT
    spider_obj.pre_run()

    # Run method
    spider_obj.run(args.gif_script, args.dbtemplate, args.working_dir)

    # Finish method to copy results
    spider_obj.finish()
