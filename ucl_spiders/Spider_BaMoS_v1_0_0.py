"""Spider_BaMoS.

Author:         Carole Sudre & Benjamin Yvernault
contact:        byvernault@ucl.ac.uk
Spider name:    BaMoS
Spider version: 1.0.0
Creation date:  2016-12-08 11:41:14.210326
Purpose:        Spider to run BaMoS, Carole script
"""

# Python packages import
import os
import sys
import shutil
from dax import spiders, SessionSpider

__author__ = "Carole Sudre & Benjamin Yvernault"
__email__ = "byvernault@ucl.ac.uk"
__purpose__ = "Spider to run BaMoS, Carole script"
__spider_name__ = "BaMoS"
__version__ = "1.0.0"
__modifications__ = """2016-12-08 11:41:14.210326 - Original write"""


BAMOS_DEFAULT = '/home/dax/Code/scripts/bash/BaMoSGenericDax.sh'
GMATRIX_DEFAULT = '/cluster/project0/SegBiASM/DataToTryBaMoS/GMatrix4_Low3.txt'
RULE_DEFAULT = '/cluster/project0/SegBiASM/DataToTryBaMoS/GenericRule_CSF.txt'
SEG_BIASM = '/home/csudre/NiftySeg_0.9.4/build_comic2/seg-apps/Seg_BiASM'
SEG_ANALYSIS = '/home/csudre/NiftySeg_0.9.4/build_comic2/seg-apps/Seg_Analysis'
TEMPLATE_SH = """
export arrayMod={array}
export Do2=1
export Do3={hast2}
export PriorsType=GIF
export NameGMatrix={gmatrix}
export RuleFileName={rule}
export MRE=250
export OW=0.01
export JC=1
export AW=0
export Opt={suffix}
export OptSP=1
export OptTA=1
export PN={session}
export PathID={dir}
export SEG_BIASM={seg_biasm}
export SEG_ANALYSIS={seg_analysis}
export REGPATH={reg_path}
export SEGPATH={seg_path}
export FSLPATH={fsl_path}
export ICBMPATH={icbm}
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
        --bamos
        --t1
        --gif
        --t2
        --flair
        --gmatrix
        --rules
        --reg_aladin
        --reg_f3d
        --seg_maths
        --seg_biasm
        --seg_analysis

    :return: argument parser object created by parse_args()
    """
    ap = spiders.get_session_argparser("BaMoS", __purpose__)
    ap.add_argument("--bamos", dest="bamos", default=BAMOS_DEFAULT,
                    help="BaMoS bash script path.")
    ap.add_argument("--t1", dest="t1", required=True,
                    help="T1 scan ID on XNAT.")
    ap.add_argument("--gif", dest="gif", required=True,
                    help="GIF proctype on XNAT.")
    ap.add_argument("--flair", dest="flair", required=True,
                    help="Flair scan ID on XNAT.")
    ap.add_argument("--t2", dest="t2", default=None,
                    help="T2 scan ID on XNAT. (optional)")
    ap.add_argument("--gmatrix", dest="gmatrix", default=GMATRIX_DEFAULT,
                    help="File defining the g-matrix.")
    ap.add_argument("--rules", dest="rulefile", default=RULE_DEFAULT,
                    help="File defining the rules.")
    ap.add_argument("--regfolder", dest="regfolder", default=None,
                    help="Folder containing the reg_ executables.")
    ap.add_argument("--segfolder", dest="segfolder", default=None,
                    help="Folder containing the seg_ executables.")
    ap.add_argument("--seg_biasm", dest="seg_biasm", default=SEG_BIASM,
                    help="Path to executable Seg_BiASM.")
    ap.add_argument("--seg_analysis", dest="seg_analysis",
                    default=SEG_ANALYSIS,
                    help="Path to executable Seg_Analysis.")
    return ap.parse_args()


class Spider_BaMoS(SessionSpider):
    """Session Spider: Spider_BaMoS.

    :param spider_path: spider file path
    :param jobdir: directory for temporary files
    :param xnat_project: project ID on XNAT
    :param xnat_subject: subject label on XNAT
    :param xnat_session: experiment label on XNAT
    :param bamos: BaMoS bash script
    :param t1: t1 scan on XNAT
    :param gif: gif proctype on XNAT
    :param t2: t2 scan on XNAT
    :param flair: flair scan on XNAT
    :param regfolder: executable for reg_*
    :param segfolder: executable for seg_*
    :param seg_biasm: executable for seg_biasm
    :param seg_analysis: executable for seg_analysis
    :param NameGMatrix: file for name g matrix
    :param RuleFileName: file for rules name
    :param xnat_host: host for XNAT if not set in environment variables
    :param xnat_user: user for XNAT if not set in environment variables
    :param xnat_pass: password for XNAT if not set in environment variables
    :param suffix: suffix to the assessor creation
    """

    def __init__(self, spider_path, jobdir,
                 xnat_project, xnat_subject, xnat_session, bamos,
                 gmatrix, rulefile, t1, gif, flair,
                 t2=None, regfolder=None, segfolder=None,
                 exe_seg_biasm=SEG_BIASM, exe_seg_analysis=SEG_ANALYSIS,
                 xnat_host=None, xnat_user=None, xnat_pass=None, suffix=""):
        """Entry point for Spider_BaMoS Class."""
        super(Spider_BaMoS,
              self).__init__(spider_path, jobdir,
                             xnat_project, xnat_subject, xnat_session,
                             xnat_host, xnat_user, xnat_pass,
                             suffix)
        self.inputs = dict()
        # XNAT extra:
        self.t1 = t1
        self.gif = gif
        self.t2 = t2
        self.flair = flair
        # args:
        self.gmatrix = gmatrix
        self.rulefile = rulefile
        # executable:
        self.bamos = bamos
        self.regfolder = regfolder
        self.segfolder = segfolder
        self.exe_seg_biasm = exe_seg_biasm
        self.exe_seg_analysis = exe_seg_analysis
        # Print version for Niftyreg
        self.check_executable(os.path.join(self.regfolder, 'reg_aladin'),
                              'reg_aladin')
        # Print version for Niftyseg
        self.check_executable('seg_LoAd', 'seg_LoAd')
        self.pdf_final = os.path.join(self.jobdir,
                                      '%sBaMoS.pdf' % self.xnat_session)

    def pre_run(self):
        """Method to download data from XNAT.

        Inputs: T1 with GIF results
                T2 (optional)
                Flair (optional)

        """
        folder = os.path.join(self.jobdir, 'inputs')
        os.makedirs(folder)
        bamos_dir = os.path.join(self.jobdir, self.xnat_session)
        os.makedirs(bamos_dir)
        os.makedirs(os.path.join(bamos_dir, 'GIF'))
        # T1:
        t1_folder = os.path.join(folder, 't1')
        os.makedirs(t1_folder)
        self.inputs['t1'] = self.download(self.t1, 'NIFTI', t1_folder)
        # GIF associated: maybe give SEG as a parameter?
        gif_label = '-x-'.join([self.xnat_project,
                                self.xnat_subject,
                                self.xnat_session,
                                self.t1,
                                self.gif])
        gif_folder = os.path.join(folder, 'gif')
        os.makedirs(gif_folder)
        self.inputs['gif'] = {}
        self.inputs['gif']['label'] = self.download(gif_label, 'LABELS',
                                                    gif_folder)
        self.inputs['gif']['tiv'] = self.download(gif_label, 'TIV',
                                                  gif_folder)
        self.inputs['gif']['prior'] = self.download(gif_label, 'PRIOR',
                                                    gif_folder)
        # Flair:
        flair_folder = os.path.join(folder, 'flair')
        os.makedirs(flair_folder)
        self.inputs['flair'] = self.download(self.flair, 'NIFTI', flair_folder)
        # Move inputs:
        shutil.move(self.inputs['t1'][0], os.path.join(
                        bamos_dir, 'T1_%s.nii.gz' % self.xnat_session))
        shutil.move(self.inputs['flair'][0], os.path.join(
                        bamos_dir, 'FLAIR_%s_init.nii.gz' % self.xnat_session))
        filelist = self.inputs['gif']['prior'] + \
            self.inputs['gif']['label'] + \
            self.inputs['gif']['tiv']
        for fpath in filelist:
            shutil.move(fpath, os.path.join(bamos_dir, 'GIF'))
        # T2
        if self.t2:
            t2_folder = os.path.join(folder, 't2')
            os.makedirs(t2_folder)
            self.inputs['t2'] = self.download(self.t2, 'NIFTI', t2_folder)
            shutil.move(self.inputs['t2'][0], os.path.join(
                    bamos_dir, 'T2_%s_init.nii.gz' % self.xnat_session))

    def run(self):
        """Method running the process for the spider on the inputs data."""
        array = '(T1 FLAIR)'
        if self.t2:
            array = '(T1 FLAIR T2)'
        args_str = TEMPLATE_SH.format(
            array=array,
            hast2=1 if self.t2 else 0,
            gmatrix=self.gmatrix,
            rule=self.rulefile,
            suffix=self.suffix,
            session=self.xnat_session,
            dir=self.jobdir,
            seg_biasm=self.exe_seg_biasm,
            seg_analysis=self.exe_seg_analysis,
            reg_path=self.regfolder,
            seg_path=self.segfolder,
            fsl_path='/share/apps/fsl-5.0.8/bin',
            icbm='/cluster/project0/SegBiASM/ICBM_Priors'
        )
        args_file = os.path.join(self.jobdir,
                                 'argument_BaMoS_%s.sh' % self.xnat_session)
        with open(args_file, "w") as f:
            f.writelines(args_str)
        cmd = "%s %s sh" % (self.bamos, args_file)
        os.system(cmd)

    def finish(self):
        """Method to copy the results in dax.RESULTS_DIR."""
        results_dict = {'PDF': self.pdf_final,
                        #
                        # ADD OTHER RESULTS YOU WANT TO SAVE
                        #
                        }
        self.upload_dict(results_dict)
        self.end()


if __name__ == '__main__':
    args = parse_args()
    # generate spider object:
    spider_obj = Spider_BaMoS(
                    spider_path=sys.argv[0],
                    jobdir=args.temp_dir,
                    xnat_project=args.proj_label,
                    xnat_subject=args.subj_label,
                    xnat_session=args.sess_label,
                    bamos=args.bamos,
                    t1=args.t1,
                    gif=args.gif,
                    t2=args.t2,
                    flair=args.flair,
                    gmatrix=args.gmatrix,
                    rulefile=args.rulefile,
                    regfolder=args.regfolder,
                    segfolder=args.segfolder,
                    exe_seg_biasm=args.seg_biasm,
                    exe_seg_analysis=args.seg_analysis,
                    xnat_host=args.host,
                    xnat_user=args.user,
                    xnat_pass=None,
                    suffix=args.suffix)
    # print some information before starting
    spider_obj.print_init(args, "Carole Sudre & Benjamin Yvernault",
                          "byvernault@ucl.ac.uk")

    # Pre-run method to download data from XNAT
    spider_obj.pre_run()

    # Run method
    spider_obj.run()

    # Finish method to copy results
    # spider_obj.finish()
