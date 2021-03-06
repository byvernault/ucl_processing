"""Spider_Verdict.

Author:         Benjamin Yvernault
contact:        b.yvernault@ucl.ac.uk
Spider name:    Verdict
Spider version: 2.0.0
Creation date:  2016-08-10 17:48:51.701453
Purpose:        Generate Verdict Map from all Verdict scans registered \
                together and merge as one nifti
"""

# Python packages import
import os
import sys
import shutil
from dax import XnatUtils, spiders, SessionSpider

__author__ = "Benjamin Yvernault"
__email__ = "b.yvernault@ucl.ac.uk"
__purpose__ = "Generate Verdict Map from all Verdict scans registered together \
and merge as one nifti"
__spider_name__ = "Verdict"
__version__ = "2.0.0"
__modifications__ = """2016-12-12 12:03 - Original write"""


PSUFFIX = 'version 2'
DEFAULT_MODEL = "VERDICTPROSTATEuDIFF1"
DEFAULT_SCHEME_FILE = "NOptimisedV_IN_v2.scheme"
DEFAULT_VERDICT_TEMPLATE = """
%% Add path from matlab folder:
addpath(genpath('{matlab_code}'));

%% Launch AMICO code:
display('Running VERDICT MAP script using AMICO');
input_path = '{input_path}';
output = '{output}';
subject = '{subject}';
filename = '{filename}';
project = '{project}';
amico = '{matlab_code}/AMICO/matlab/';
camino = '{camino}';
spams = '{spams}';
scheme = '{scheme_filename}';
model = '{model}';
launch_AMICO_for_INNOVATE(input_path,subject,filename,output,project,amico,\
camino,spams,scheme,model);

%% Generate DICOM:
nb_acq = {acq};
nii_folder = [output '/AMICO/' model];
dicom_file = '{dicom_file}';
matlab_nii_code = '{matlab_code}/ext';
out_dcm = '{out_dcm}';
suffix = '{suffix}';

display('Converting NIFTI to DICOM in color.');
nii2dicomRGB(nii_folder, dicom_file, out_dcm, matlab_nii_code, nb_acq, suffix)
"""
PDF_TEMPLATE = """
%% Add path from matlab folder:
addpath(genpath('{matlab_code}'));
%% Generate PDF for VERDICT MAP
% Maps Name:
model = '{model}';
nii_folder = ['{output}' '/AMICO/' model];
subject = '{subject}';
output_pdf = '{out_pdf}';
nb_acq = {acq};
maps_name = {{'fIC','R','cellularity','fEES','fVASC','FobjCamino'}};

% Open an image to see the number of slices
f_file = fullfile(nii_folder, 'FIT_FobjCamino.nii');
aux = load_untouch_nii(f_file);
nb_slices = aux.hdr.dime.dim(4);
display('Generating the PDF for VERDICT maps.');
for i=1:nb_slices
    filename = ['DisplayMapsVerdictSlice' num2str(i, '%02d');];
    jpg_path = fullfile(output_pdf, [filename '.pdf']);
    plot_oneslice_selectedmaps(nii_folder,subject,i,maps_name,1,1,1,jpg_path,nb_acq);
end
"""
DICOM_SCAN_TYPE = ['WIP b3000_90 SENSE', 'SWITCH DB TO YES b3000_80',
                   'b3000_80']


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
        --proctype: proctype where inputs are
        --mc: matlab code to launch verdict map
        --amico : path to AMICO folder
        --camino : path to Camino
        --spams : path to spams-matlab

    :return: argument parser object created by parse_args()
    """
    ap = spiders.get_session_argparser("Verdict", __purpose__)
    ap.add_argument("--proctype", dest="proctype", default=None, required=True,
                    help="Assessor type containing the registered nifti.")
    ap.add_argument("--nbAcq", dest="nb_acquisition", default=1,
                    required=True,
                    help="Number of Acquisition of VERDICT scans (1 or 2).")
    ap.add_argument("--mc", dest="matlab_code", default=None, required=True,
                    help="Matlab code folder where is \
launch_AMICO_for_INNOVATE and the AMICO folder.")
    ap.add_argument("--camino", dest="camino", default=None, required=True,
                    help="Path to Camino folder.")
    ap.add_argument("--spams", dest="spams", default=None, required=True,
                    help="Path to spams-matlab folder.")
    ap.add_argument("--scheme", dest="scheme_filename", default=None,
                    help="Path to schemeFilename \
(NOptimisedV_IN.scheme).")
    ap.add_argument("--model", dest="model", default=DEFAULT_MODEL,
                    help="Model type for AMICO: VERDICTPROSTATEuDIFF1.")
    return ap.parse_args()


class Spider_Verdict(SessionSpider):
    """Session Spider: Spider_Verdict.

    :param spider_path: spider file path
    :param jobdir: directory for temporary files
    :param xnat_project: project ID on XNAT
    :param xnat_subject: subject label on XNAT
    :param xnat_session: experiment label on XNAT
    :param matlab_path: path to matlab code
    :param xnat_host: host for XNAT if not set in environment variables
    :param xnat_user: user for XNAT if not set in environment variables
    :param xnat_pass: password for XNAT if not set in environment variables
    :param suffix: suffix to the assessor creation
    """

    def __init__(self, spider_path, jobdir, xnat_project, xnat_subject,
                 xnat_session, proctype, nb_acquisition,
                 matlab_code, camino, spams,
                 scheme_filename=None, model=DEFAULT_MODEL,
                 xnat_host=None, xnat_user=None, xnat_pass=None,
                 suffix=""):
        """Entry point for Spider_Verdict Class."""
        super(Spider_Verdict,
              self).__init__(spider_path, jobdir,
                             xnat_project, xnat_subject, xnat_session,
                             xnat_host, xnat_user, xnat_pass,
                             suffix)
        self.proctype = proctype
        self.nb_acquisition = int(nb_acquisition)

        self.inputs = dict()
        self.pdf_final = os.path.join(self.jobdir,
                                      '%s_VERDICT_report.pdf' % xnat_session)

        self.matlab_code = matlab_code
        self.spams = spams
        self.camino = camino
        if not scheme_filename:
            self.scheme_filename = os.path.join(matlab_code,
                                                DEFAULT_SCHEME_FILE)
        else:
            self.scheme_filename = scheme_filename
        self.model = model

    def pre_run(self):
        """Method to download data from XNAT.

        :param argument_parse: argument parser object return by parse_args()
        """
        resource = 'ACQ'  # resource to download from the scan on XNAT
        folder = os.path.join(self.jobdir, 'inputs')
        os.makedirs(folder)
        assessor_label = '-x-'.join([self.xnat_project,
                                     self.xnat_subject,
                                     self.xnat_session,
                                     self.proctype])

        xnat = XnatUtils.get_interface()
        a = xnat.select('/projects/%s/subjects/%s/experiments/%s/assessors/%s'
                        % (self.xnat_project, self.xnat_subject,
                           self.xnat_session, assessor_label))

        for nb_acq in range(1, self.nb_acquisition+1):
            res_name = '%s%d' % (resource, nb_acq)
            self.time_writer('Download resource: %s' % res_name)
            file_gzip = XnatUtils.download_file_from_obj(
                                  folder, a.out_resource(res_name))
            self.time_writer('Unzip file: %s' % file_gzip)
            XnatUtils.gunzip_file(file_gzip)
            self.inputs[nb_acq] = file_gzip[:-3]

        sc = XnatUtils.get_good_scans(a.parent(), DICOM_SCAN_TYPE)[0]
        dcm_file = XnatUtils.download_file_from_obj(
                              folder, sc.resource('DICOM'))
        self.inputs['dcm'] = dcm_file
        xnat.disconnect()

    def run(self):
        """Method running the process for the spider on the inputs data."""
        output_folder = XnatUtils.makedir(os.path.join(self.jobdir, 'outputs'))
        dcm_folder = XnatUtils.makedir(os.path.join(output_folder, 'OsiriX'))
        pdfs_dir = XnatUtils.makedir(os.path.join(output_folder, 'pdfs'))
        fpages = list()
        # Load dicom headers
        if not os.path.isfile(self.inputs['dcm']):
            err = "DICOM File %s not found."
            raise Exception(err % self.inputs['dcm'])

        for nb_acq in range(1, self.nb_acquisition+1):
            folder = os.path.join(output_folder, str(nb_acq))
            os.makedirs(folder)
            mat_lines = DEFAULT_VERDICT_TEMPLATE.format(
                    matlab_code=self.matlab_code,
                    input_path=os.path.dirname(self.inputs[nb_acq]),
                    subject=self.xnat_subject,
                    filename=os.path.basename(self.inputs[nb_acq]),
                    output=folder,
                    project=self.xnat_project,
                    camino=self.camino,
                    spams=self.spams,
                    scheme_filename=self.scheme_filename,
                    model=self.model,
                    dicom_file=self.inputs['dcm'],
                    out_dcm=dcm_folder,
                    acq=str(nb_acq),
                    suffix=PSUFFIX,
                    )
            matlab_script = os.path.join(output_folder,
                                         'run_verdict_map%d.m' % nb_acq)
            with open(matlab_script, "w") as f:
                f.writelines(mat_lines)
            self.run_matlab(matlab_script, verbose=True)

            mat_lines = PDF_TEMPLATE.format(
                    matlab_code=self.matlab_code,
                    subject=self.xnat_subject,
                    output=folder,
                    model=self.model,
                    out_pdf=pdfs_dir,
                    acq=str(nb_acq),
                    )
            matlab_script = os.path.join(output_folder,
                                         'run_pdf_%d.m' % nb_acq)
            with open(matlab_script, "w") as f:
                f.writelines(mat_lines)
            XnatUtils.run_matlab(matlab_script, verbose=True)
            pdf_pages = XnatUtils.find_files(pdfs_dir, '.pdf')
            # Merge all pdfs into one:
            pdf_page = os.path.join(output_folder, str(nb_acq),
                                    'VerdictMapAcq%d.pdf' % nb_acq)
            self.merge_pdf_pages(pdf_pages, pdf_page)
            fpages.append(pdf_page)

        # Merge PDFs:
        if len(fpages) > 1:
            self.merge_pdf_pages(fpages, self.pdf_final)
        else:
            shutil.move(fpages[0], self.pdf_final)

        # Zip the DICOMs output:
        initdir = os.getcwd()
        # Zip all the files in the directory
        zip_name = os.path.join(self.jobdir, 'outputs', 'OsiriX', 'osirix.zip')
        os.chdir(os.path.join(self.jobdir, 'outputs', 'OsiriX'))
        os.system('zip -r %s * > /dev/null' % zip_name)
        # return to the initial directory:
        os.chdir(initdir)

        # Gzip nii:
        XnatUtils.gzip_nii(os.path.join(output_folder, '1', 'AMICO',
                           self.model))
        if self.nb_acquisition == 2:
            XnatUtils.gzip_nii(os.path.join(output_folder, '2', 'AMICO',
                                            self.model))

    def finish(self):
        """Method to copy the results in dax.RESULTS_DIR."""
        results_dict = {'PDF': self.pdf_final,
                        'OsiriX': os.path.join(self.jobdir, 'outputs',
                                               'OsiriX', 'osirix.zip')}
        for nb_acq in range(1, self.nb_acquisition+1):
            res = 'RMAPS%d' % nb_acq
            results_dict[res] = os.path.join(self.jobdir, 'outputs',
                                             str(nb_acq), 'AMICO', self.model)
        self.upload_dict(results_dict)
        self.end()

    def run_matlab(self, matlab_script, verbose=False):
        """Call MATLAB with -nodesktop -nosplash and -singlecompthread.

        :param matlab_script: Full path to the .m file to run
        :param verbose: True to print all MATLAB output to terminal, False to
         suppress.
        :return: None
        """
        self.time_writer("Matlab script: %s running ..." % matlab_script)
        cmd = "matlab -nodisplay -nodesktop -nojvm -nosplash -singleCompThread \
    < %s" % matlab_script
        if not verbose:
            matlabdir = os.path.dirname(matlab_script)
            prefix = os.path.basename(matlab_script).split('.')[0]
            cmd = cmd+' > '+os.path.join(matlabdir, prefix+'_outlog.log')
        os.system(cmd)
        self.time_writer("Matlab script: %s done" % matlab_script)


if __name__ == '__main__':
    args = parse_args()
    # generate spider object:
    spider_obj = Spider_Verdict(spider_path=sys.argv[0],
                                jobdir=args.temp_dir,
                                xnat_project=args.proj_label,
                                xnat_subject=args.subj_label,
                                xnat_session=args.sess_label,
                                proctype=args.proctype,
                                nb_acquisition=args.nb_acquisition,
                                matlab_code=args.matlab_code,
                                camino=args.camino,
                                spams=args.spams,
                                scheme_filename=args.scheme_filename,
                                model=args.model,
                                xnat_host=args.host,
                                xnat_user=args.user,
                                xnat_pass=None,
                                suffix=args.suffix)

    # print some information before starting
    spider_obj.print_init(args, "Benjamin Yvernault", "b.yvernault@ucl.ac.uk")

    # Pre-run method to download data from XNAT
    spider_obj.pre_run()

    # Run method
    spider_obj.run()

    # Finish method to copy results
    if not args.skipfinish:
        spider_obj.finish()
