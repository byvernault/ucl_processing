from dax import XnatUtils, spiders
import os
import shutil


# Inputs:
LABEL = '${assessor_label}'
JOB_DIR = '${temp_dir}'
ACQ1_FILE = '${acquisition1}'
ACQ2_FILE = '${acquisition2}'

# Variables
DEFAULT_MODEL = "VerdictProstate_Rmaps"
DEFAULT_SCHEME_FILE = "NOptimisedV_IN.scheme"
DEFAULT_VERDICT_TEMPLATE = """
%% Add path from matlab folder:
addpath(genpath('${matlab_src}'));

%% Launch AMICO code:
display('Running VERDICT MAP script using AMICO');
input_path = '{input_path}';
output = '{output}';
subject = '{subject}';
filename = '{filename}';
project = '{project}';
amico = '${matlab_src}/AMICO/matlab/';
camino = '${camino}';
spams = '${spams}';
scheme = '${scheme_filename}';
model = '${model}';
launch_AMICO_for_INNOVATE(input_path,subject,filename,output,project,amico,\
camino,spams,scheme,model);

%% Generate DICOM:
nb_acq = {acq};
nii_folder = [output '/AMICO/' model];
dicom_file = '${dcm_file}';
matlab_nii_code = '${matlab_src}/ext';
out_dcm = '{out_dcm}';
suffix = 'version 1';

display('Converting NIFTI to DICOM in color.');
nii2dicomRGB(nii_folder, dicom_file, out_dcm, matlab_nii_code, nb_acq, suffix)
"""

PDF_TEMPLATE = """
%% Add path from matlab folder:
addpath(genpath('${matlab_src}'));
%% Generate PDF for VERDICT MAP
% Maps Name:
model = '${model}';
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


def main():
    """ Main function."""
    if len(JOB_DIR) > 0 and not os.path.exists(JOB_DIR):
        os.makedirs(JOB_DIR)

    dcm_folder = XnatUtils.makedir(os.path.join(JOB_DIR, 'OsiriX'))
    pdfs_dir = XnatUtils.makedir(os.path.join(JOB_DIR, 'pdfs'))
    fpages = list()

    project, subject, session = LABEL.split('-x-')[:3]
    pdf_final = os.path.join(JOB_DIR, '%s_VERDICT_report.pdf' % session)

    for nb_acq, acq_path in enumerate([ACQ1_FILE, ACQ2_FILE]):
        nb_acq = nb_acq + 1
        if not os.path.exists(acq_path):
            continue
        else:
            folder = os.path.join(JOB_DIR, str(nb_acq))
            os.makedirs(folder)

            # Unzip niftis:
            XnatUtils.gunzip_file(acq_path)
            acq_path = acq_path[:-3]

            # Command
            mat_lines = DEFAULT_VERDICT_TEMPLATE.format(
                input_path=os.path.dirname(acq_path),
                filename=os.path.basename(acq_path),
                output=folder,
                out_dcm=dcm_folder,
                acq=str(nb_acq),
                project=project,
                subject=subject,
            )
            matlab_script = os.path.join(JOB_DIR,
                                         'run_verdict_map%d.m' % nb_acq)
            with open(matlab_script, "w") as f:
                f.writelines(mat_lines)
            cmd = "matlab -nodisplay -nodesktop -nosplash -singleCompThread < "
            os.system(cmd + matlab_script)

            mat_lines = PDF_TEMPLATE.format(
                output=folder,
                out_pdf=pdfs_dir,
                acq=str(nb_acq),
                subject=subject,
            )
            matlab_script = os.path.join(JOB_DIR, 'run_pdf_%d.m' % nb_acq)
            with open(matlab_script, "w") as f:
                f.writelines(mat_lines)
            XnatUtils.run_matlab(matlab_script, verbose=True)

            pdf_pages = XnatUtils.find_files(pdfs_dir, '.pdf')
            # Merge all pdfs into one:
            pdf_page = os.path.join(JOB_DIR, str(nb_acq),
                                    'VerdictMapAcq%d.pdf' % nb_acq)
            spiders.merge_pdfs(pdf_pages, pdf_page)
            fpages.append(pdf_page)

    # Merge PDFs:
    if len(fpages) > 1:
        spiders.merge_pdfs(fpages, pdf_final)
    else:
        shutil.move(fpages[0], pdf_final)

    # Zip the DICOMs output:
    initdir = os.getcwd()
    # Zip all the files in the directory
    zip_name = os.path.join(JOB_DIR, 'OsiriX', 'osirix.zip')
    os.chdir(os.path.join(JOB_DIR, 'OsiriX'))
    os.system('zip -r %s * > /dev/null' % zip_name)
    # return to the initial directory:
    os.chdir(initdir)

    # Gzip nii:
    XnatUtils.gzip_nii(os.path.join(JOB_DIR, '1', 'AMICO', '${model}'))
    if ACQ2_FILE and os.path.exists(ACQ2_FILE):
        XnatUtils.gzip_nii(os.path.join(JOB_DIR, '2', 'AMICO', '${model}'))


if __name__ == '__main__':
    main()
