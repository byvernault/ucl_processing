"""Spider to perform Diffusion Model Fitting.

Author:         Benjamin Yvernault
contact:        b.yvernault@ucl.ac.uk
Spider name:    Diffusion_Model_Fitting
Spider version: 1.0.0
Creation date:  2016-05-10 13:36:26.321786
Purpose:        Perform Diffusion Model Fitting with pre-processing steps.
"""

# Python packages import
import os
import sys
import subprocess as sb
from dax import spiders, SessionSpider

__author__ = "Benjamin Yvernault"
__email__ = "byvernault@gmail.com"
__purpose__ = "Perform Diffusion Model Fitting with pre-processing steps."
__spider_name__ = "Diffusion_Model_Fitting"
__version__ = "1.0.0"
__modifications__ = """2016-05-10 13:36:26.321786 - Original write"""

DTI_CMD = """python {exe_path} \
-i {images} \
-a {bvals} \
-e {bvecs} \
-t {gif_t1} \
--t1_mask {gif_brain} \
-o {OUTPUT_FOLDER} \
--no_qsub \
--n_procs 1 \
--openmp_core {number_core} \
--working_dir '{working_dir}' \
{args}
"""
DEFAULT_ARGS = "--rot 34.56 --etd 2.46 --ped -y"


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
        --dti : dti scan IDs on XNAT
        --gif  : Label of GIF assessor on XNAT
        --fmagni : Field Map Magnitude image file to be associated with
                   the DWIs
        --fphase : Field Map Phase image file to be associated with the
                   DWIs
        --number_core : number of core used by the process
        --exe : path to the executable perform_dti_processing.py
        --dtiargs : other arguments for perform_dti_processing as a string
                    (rot/etd/ped). Eg: --rot 34.56 --etd 2.46 --ped -y
        --working_dir : working directory for temp files
    :return: argument parser object created by parse_args()
    """
    ap = spiders.get_session_argparser("Diffusion_Model_Fitting", __purpose__)
    ap.add_argument("--dti", dest="dtis", required=True, metavar='sep ","',
                    help="ID(s) on XNAT for dti scan.")
    ap.add_argument("--gif", dest="gif", required=True,
                    help="Label of GIF assessor on XNAT.")
    ap.add_argument("--fmagni", dest="fmagni", Default=None,
                    help="Field Map Magnitude image file to be associated with\
                    the DWIs.")
    ap.add_argument("--fphase", dest="fphase", Default=None,
                    help="Field Map Phase image file to be associated with the\
                    DWIs")
    ap.add_argument("--exe", dest="dti_exe", required=True,
                    help="Path to perform_dti_processing.py.")
    ap.add_argument("--exeargs", dest="dti_args", help="other arguments for \
perform_dti_processing as a string (rot/etd/ped). \
Default: --rot 34.56 --etd 2.46 --ped -y",
                    default="--rot 34.56 --etd 2.46 --ped -y")
    ap.add_argument("--openmp_core", dest="openmp_core",
                    help="Number of core used.", default=1)
    ap.add_argument("--working_dir", dest="working_dir", default="",
                    help="working directory for temp files. Default: output")
    return ap.parse_args()


class Spider_Diffusion_Model_Fitting(SessionSpider):
    """Spider to perform Diffusion Model Fitting with pre-processing steps.

    :param spider_path: spider file path
    :param jobdir: directory for temporary files
    :param number_core: number of core used by reg_aladin
    :param xnat_project: project ID on XNAT
    :param xnat_subject: subject label on XNAT
    :param xnat_session: experiment label on XNAT
    :param xnat_host: host for XNAT if not set in environment variables
    :param xnat_user: user for XNAT if not set in environment variables
    :param xnat_pass: password for XNAT if not set in environment variables
    :param suffix: suffix to the assessor creation
    """

    def __init__(self, spider_path, jobdir,
                 dti_exe, dtiargs, dti, gif, fmagni, fphase,
                 xnat_project, xnat_subject, xnat_session,
                 xnat_host=None, xnat_user=None, xnat_pass=None,
                 number_core=1, suffix=""):
        """Init Method."""
        super(Spider_Diffusion_Model_Fitting, self).__init__(spider_path,
                                                             jobdir,
                                                             xnat_project,
                                                             xnat_subject,
                                                             xnat_session,
                                                             xnat_host,
                                                             xnat_user,
                                                             xnat_pass,
                                                             suffix)
        self.inputs = {'dti': list(),
                       't1': '',
                       'gif': ''}
        self.dti_exe = dti_exe
        self.dtiargs = dtiargs
        self.dti = dti
        self.gif = gif
        self.fmagni = fmagni
        self.fphase = fphase
        self.number_core = number_core
        # Print version for Niftyreg  (NOT CHECKED EXCEPT REG_ALADIN)
        # List of the binaries necessary for this pipeline:
        # * FSL: fslmaths, fslmerge, fslsplit
        # * niftyreg: reg_aladin, reg_transform, reg_f3d, reg_resample
        # * niftyseg: seg_maths
        # * niftyfit: fit_dwi, dwi_tool
        # * susceptibility: pm_scale, gen_fm, gen_pm, phase_unwrap
        exe_call = sb.Popen(['reg_aladin', '--version'],
                            stdout=sb.PIPE,
                            stderr=sb.PIPE)
        exe_version, _ = exe_call.communicate()
        self.time_writer('Niftyreg version: %s' % (exe_version.strip()))
        self.time_writer('GIF version: 3a76a255ab')

    def pre_run(self, argument_parse):
        """Method to download data from XNAT.

        :param argument_parse: argument parser object return by parse_args()
        """
        folder = os.path.join(self.jobdir, 'inputs')
        os.makedirs(folder)
        # dti(s) with bval and bvec
        for scan_id in self.dti.split(','):
            dti_dir = os.path.join(folder, scan_id)
            os.makedirs(dti_dir)
            image = self.download(scan_id, 'NIFTI', dti_dir)
            if not image:
                raise Exception("No NIFTI downloaded from XNAT for %s." %
                                (scan_id))
            bval = self.download(scan_id, 'BVAL', dti_dir)
            if not bval:
                raise Exception("No BVAL downloaded from XNAT for %s." %
                                (scan_id))
            bvec = self.download(scan_id, 'BVEC', dti_dir)
            if not bvec:
                raise Exception("No BVEC downloaded from XNAT for %s." %
                                (scan_id))
            self.inputs['dti'].append({'image': image,
                                       'bval': bval,
                                       'bvec': bvec})
        # T1 scan
        self.inputs['t1'] = self.download(self.gif.split('-x-')[-2],
                                          'NIFTI', folder)
        if not self.inputs['t1']:
            raise Exception("No NIFTI downloaded from XNAT for %s." %
                            (scan_id))
        # GIF results:
        gif = '-x-'.join([self.xnat_project,
                          self.xnat_subject,
                          self.xnat_session,
                          self.t1,
                          'GIF_Parcellation_v1'])
        self.inputs['gif'] = self.download(gif, 'BRAIN', folder)
        if not self.inputs['gif']:
            raise Exception("No BRAIN resource downloaded from XNAT for GIF.")

    def run(self, exe_path, working_dir):
        """Method running the process for the spider on the inputs data."""
        if len(working_dir) > 0 and not os.path.exists(working_dir):
            os.makedirs(working_dir)

        if not exe_path.endswith('perform_gif_propagation.py'):
            exe_path = os.path.join(exe_path, "perform_gif_propagation.py")

        if not os.path.exists(exe_path):
            raise Exception("Python script: %s not found" % (exe_path))
        else:
            if not os.path.exists(os.path.join(self.jobdir, 'outputs')):
                os.makedirs(os.path.join(self.jobdir, 'outputs'))
            dtis = ' '.join([d('image') for d in self.inputs['dti']])
            bvals = ' '.join([d('bval') for d in self.inputs['dti']])
            bvecs = ' '.join([d('bvec') for d in self.inputs['dti']])
            cmd = DTI_CMD.format(exe_path=exe_path,
                                 images=dtis,
                                 bvals=bvals,
                                 bvecs=bvecs,
                                 gif_t1=self.inputs['t1'],
                                 gif_brain=self.inputs['gif'],
                                 output=os.path.join(self.jobdir, 'outputs'),
                                 number_core=self.number_core,
                                 working_dir=working_dir,
                                 args=self.dtiargs)
            if self.fmagni:
                cmd = "%s -m %s" % (cmd, self.fmagni)
            if self.fphase:
                cmd = "%s -p %s" % (cmd, self.fphase)
            self.run_system_cmd(cmd)
            # self.make_pdf()

    def finish(self):
        """Method to copy the results in the dax.RESULTS_DIR folder."""
        self.time_writer('Results saved in dax.RESULTS_DIR')
        results_dict = {'PDF': ''}
        self.upload_dict(results_dict)
        self.end()

if __name__ == '__main__':
    args = parse_args()
    # generate spider object:
    spider_obj = Spider_Diffusion_Model_Fitting(spider_path=sys.argv[0],
                                                jobdir=args.temp_dir,
                                                dti_exe=args.dti_exe,
                                                dtiargs=args.dtiargs,
                                                dti=args.dti,
                                                gif=args.gif,
                                                fmagni=args.fmagni,
                                                fphase=args.fphase,
                                                xnat_project=args.proj_label,
                                                xnat_subject=args.subj_label,
                                                xnat_session=args.sess_label,
                                                xnat_scan=args.scan_label,
                                                xnat_host=args.host,
                                                xnat_user=args.user,
                                                xnat_pass=None,
                                                suffix=args.suffix)
    # print some information before starting
    spider_obj.print_init(args, "Benjamin Yvernault", "byvernault@gmail.com")

    # Pre-run method to download data from XNAT
    spider_obj.pre_run(args)

    # Run method
    spider_obj.run(args.exe_path, args.working_dir)

    # Finish method to copy results
    spider_obj.finish()
