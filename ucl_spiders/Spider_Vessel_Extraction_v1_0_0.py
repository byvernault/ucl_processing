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
import subprocess as sb
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
VESSEL_CMD = """{exe_path} -i {input} -o {output} --mod 0 --aone 0.5 --atwo 2\
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
    :param suffix: suffix to the assessor creation
    """

    def __init__(self, spider_path, jobdir, xnat_project, xnat_subject,
                 xnat_session, xnat_scan, xnat_host=None, xnat_user=None,
                 xnat_pass=None, suffix=""):
        """Init function for the Spider class."""
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
        self.pixel_size = DEFAULT_PIXEL_SIZE
        # Print version for Niftyreg - GIFi
        pversion = sb.Popen([ARGS.vessel_extractor, '--version'],
                            stdout=sb.PIPE,
                            stderr=sb.PIPE)
        nve_version, _ = pversion.communicate()
        self.time_writer('niftkVesselExtractor version: %s' %
                         (nve_version.strip()))

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

    def run(self, vessel_extractor):
        """Method running the process for the spider on the inputs data.

        :return: None
        """
        if vessel_extractor.endswith('niftkVesselExtractor'):
            exe_path = vessel_extractor
        elif os.path.isdir(vessel_extractor):
            exe_path = os.path.join(vessel_extractor, "niftkVesselExtractor")

        if not os.path.exists(exe_path):
            raise Exception("Executable '%s' not found" % (exe_path))
        else:
            if not os.path.exists(os.path.join(self.jobdir, 'outputs')):
                os.makedirs(os.path.join(self.jobdir, 'outputs'))
            output_name = ("%s_vessel_extracted.nii.gz" %
                           (os.path.basename(self.inputs[0]).split('.')[0]))
            self.output = os.path.join(self.jobdir, 'outputs', output_name)
            cmd = VESSEL_CMD.format(exe_path=exe_path,
                                    input=self.inputs[0],
                                    output=self.output,
                                    pixel_size=self.pixel_size)
            self.run_system_cmd(cmd)
            # self.make_pdf()

    def finish(self):
        """Method to copy the results in dax.RESULTS_DIR.

        :return: None
        """
        results_dict = {
                        # 'PDF': pdfpath,
                        # 'OUTPUT': self.output
                        }
        self.upload_dict(results_dict)
        self.end()


if __name__ == '__main__':
    ARGS = parse_args()
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
                                          suffix=ARGS.suffix)
    # print some information before starting
    spider_obj.print_init(ARGS, "Benjamin Yvernault", "b.yvernault@ucl.ac.uk")

    # Pre-run method to download data from XNAT
    spider_obj.pre_run()

    # Run method
    spider_obj.run(ARGS.vessel_extractor)

    # Finish method to copy results
    spider_obj.finish()
