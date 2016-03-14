""" Module to automatically archive sessions on XNAT following REDCap data"""
import os, re, time, logging

import redcap

from dax import XnatUtils, SessionModule
from VUIIS_path_settings import API_URL

LOGGER = logging.getLogger('dax')

DEFAULT_MODULE_NAME = 'Auto Archive'
DEFAULT_TMP_PATH = os.path.join('/tmp', DEFAULT_MODULE_NAME)
DEFAULT_TEXT_REPORT = 'ERROR in module ' + DEFAULT_MODULE_NAME + ':\n'
DEFAULT_API_KEY = 'API_KEY_Module_Auto_Archive'
DEFAULT_RC_SUBJ_FIELD = 'subj_id'
DEFAULT_RC_SESS_FIELD = 'sess_id'
DEFAULT_RC_PROJ_FIELD = 'project_name'
DEFAULT_RC_DATE_FIELD = 'sess_date'
DEFAULT_RC_VUIIS_FIELD = 'scan_id'
RESET_SLEEP_SECS = 30

PROJ_DICT = {'Cognitive Complaints': 'NewhouseCC',
             'Menstrual Math' : 'NewhouseMM',
             'Procedural Learning' : 'NewhousePL',
             'MDD hx' : 'NewhouseMDDHx',
             'Breast Cancer Study' : 'NewhouseBC'}

class Module_Auto_Archive(SessionModule):
    """ Module to automatically archive sessions on XNAT """
    def __init__(self,
                 mod_name=DEFAULT_MODULE_NAME,
                 directory=DEFAULT_TMP_PATH,
                 email=None,
                 text_report=DEFAULT_TEXT_REPORT,
                 api_url=API_URL,
                 api_key=DEFAULT_API_KEY,
                 rc_proj_field=DEFAULT_RC_PROJ_FIELD,
                 rc_subj_field=DEFAULT_RC_SUBJ_FIELD,
                 rc_sess_field=DEFAULT_RC_SESS_FIELD,
                 rc_date_field=DEFAULT_RC_DATE_FIELD,
                 rc_vuiis_field=DEFAULT_RC_VUIIS_FIELD,
                 rc_subj_event=None,
                 rc_sess_event=None,
                 proj_map='0',
                 pre_proj_filter=None,
                 arc_proj_filter=None):
        """ init function overridden from base-class """
        super(Module_Auto_Archive, self).__init__(mod_name, directory, email, \
                                                 text_report=text_report)

        self.api_url = api_url

        if api_key in os.environ:
            self.api_key = os.environ[api_key]
        else:
            self.api_key = api_key

        self.rc_proj_field = rc_proj_field
        self.rc_subj_field = rc_subj_field
        self.rc_sess_field = rc_sess_field
        self.rc_date_field = rc_date_field
        self.rc_vuiis_field = rc_vuiis_field

        self.rc_fields = [self.rc_proj_field, self.rc_subj_field, self.rc_sess_field, \
                          self.rc_date_field, self.rc_vuiis_field]
        self.rc_subj_event = rc_subj_event
        self.rc_sess_event = rc_sess_event
        self.proj_map = proj_map
        self.pre_proj_filter = pre_proj_filter
        self.arc_proj_filter = arc_proj_filter
        self.rc_list = ()
        self.pa_list = ()
        self.xnat = None

    def change_prearchive_project(self, pre_uri, new_project, reset_first=False):
        """ change the project assign to a session in the prearchive """

        if reset_first:
            # Force a refresh of files and wait
            LOGGER.debug('resetting session in Archive and waiting...')
            result = self.reset_prearchive_session(pre_uri)
            LOGGER.debug('RESULT=' + result)
            time.sleep(RESET_SLEEP_SECS)

        post_body = """src={src}&newProject={proj}&async=false""".format(src=pre_uri,
                                                                         proj=new_project)
        request_uri = '/data/services/prearchive/move'
        return self.xnat._exec(request_uri, 'POST', post_body,
                               {'content-type':'application/x-www-form-urlencoded'})

    def reset_prearchive_session(self, pre_uri):
        """ reset session in the prearchive """
        post_body = "action=build"
        request_uri = '/data' + pre_uri
        return self.xnat._exec(request_uri, 'POST', post_body,
                               {'content-type':'application/x-www-form-urlencoded'})

    def archive_prearchive_session(self, src, proj, subj, exp, reset_first=False):
        """ archive the session from the prearchive """

        if reset_first:
            # Force a refresh of files and wait
            LOGGER.info('resetting session in Archive and waiting...')
            result = self.reset_prearchive_session(src)
            LOGGER.debug('RESULT=' + result)
            time.sleep(RESET_SLEEP_SECS)

        post_body = """src={src}&project={proj}&subject={subj}&session={sess}""".format(src=src,
                                                                                        proj=proj,
                                                                                        subj=subj,
                                                                                        sess=exp)
        request_uri = '/data/services/archive'
        return self.xnat._exec(request_uri, 'POST', post_body,
                               {'content-type':'application/x-www-form-urlencoded'})

    def load_redcap(self):
        """ load redcap project """
        rc = redcap.Project(self.api_url, self.api_key)
        rc_sess_fields = self.rc_fields + [rc.def_field]

        if self.rc_sess_event == None:
            sess_event_list = None
        else:
            sess_event_list = [self.rc_sess_event]

        if self.rc_subj_event == None:
            subj_event_list = None
        else:
            subj_event_list = [self.rc_subj_event]

        # Load sess events
        self.rc_list = rc.export_records(fields=rc_sess_fields,
                                         raw_or_label='label',
                                         events=sess_event_list)



        # Only keep records with a vuiis number
        self.rc_list = [r for r in self.rc_list if r[self.rc_vuiis_field]]

        if self.rc_subj_event != self.rc_sess_event:
            # Load events with subj
            rc_subj_fields = [rc.def_field, self.rc_subj_field]
            rc_subj_list = rc.export_records(fields=rc_subj_fields,
                                             raw_or_label='label',
                                             events=subj_event_list)

            id2subj = dict([(x[rc.def_field], x[self.rc_subj_field]) for x in rc_subj_list])

            # Add subj id to sess events
            for r in self.rc_list:
                subj_id = r[rc.def_field]
                r[self.rc_subj_field] = id2subj[subj_id]

    def load_prearchive(self):
        """ retrieve sessions from the prearchive """
        if self.xnat != None:
            self.pa_list = self.xnat._get_json('/data/prearchive')
        else:
            LOGGER.error('cannot load prearchive, no XNAT connection.')

    def check_projects(self):
        """ check the projects between prearchive and redcap """

         # Filter the list
        if self.pre_proj_filter:
            self.pa_list = [p for p in self.pa_list if p['project'] == self.pre_proj_filter]

        # Check each session in prearchive against redcap
        for p in self.pa_list:
            pre_proj = p['project']
            pre_date = p['scan_date'][:10]  # get only the date without time
            pre_num = p['name']
            pre_uri = p['url']
            pre_stat = p['status']

            # Check prearchive status, must be Ready to archive
            if pre_stat != 'READY':
                LOGGER.info('Prearchive session not ready, status=' + pre_stat)
                continue

            # Look for a matching scan in RedCap
            rc_scan = None
            for r in self.rc_list:
                vuiis_label = str(r[self.rc_vuiis_field])

                # Remove PI prefix if found
                if re.match("[A-Za-z]+_\d+", vuiis_label):
                    vuiis_label = vuiis_label.split('_')[1]

                if vuiis_label != pre_num:
                    continue
                elif r[self.rc_date_field] != pre_date:
                    msg = 'date mismatch:{0}, REDCap={1}, XNAT={2}'.format(pre_num,
                                                                           r[self.rc_date_field],
                                                                           pre_date)
                    self.report_error(msg)
                    continue
                else:
                    rc_scan = r
                    break

            if not rc_scan:
                msg = 'no match found for Prearchive Session:' + pre_num
                self.report_error(msg)
                continue

            LOGGER.debug('found match for Prearchive Session:' + pre_num)
            for f in self.rc_fields:
                LOGGER.debug(f + '=' + str(rc_scan[f]))

            proj_label = str(rc_scan[self.rc_proj_field])
            if self.proj_map == '1':
                proj_label = PROJ_DICT[proj_label]

            if pre_proj == proj_label:  # Skip if already moved to correct project
                LOGGER.debug('project already correct')
                continue

            # Check for project on XNAT, must exist
            xnat_proj = self.xnat.select.project(proj_label)
            if not xnat_proj.exists():
                msg = 'project does not exist in XNAT:' + proj_label
                self.report_error(msg)
                continue

            # Change project
            LOGGER.info('change Project for Session:{0} to {1}...'.format(pre_uri, proj_label))
            result = self.change_prearchive_project(pre_uri, proj_label, reset_first=True)
            LOGGER.debug('RESULT=' + result)

    def do_archiving(self):
        """ archive the sessions """

        # Filter the list
        if self.arc_proj_filter:
            self.pa_list = [p for p in self.pa_list if p['project'] == self.arc_proj_filter]

        LOGGER.info('checking Prearchive Scans:n={0}'.format(len(self.pa_list)))
        for p in self.pa_list:
            pre_proj = p['project']
            pre_date = p['scan_date'][:10]  # get only the date without time
            pre_num = p['name']
            pre_uri = p['url']
            pre_stat = p['status']

            LOGGER.debug('checking Prearchive Session:project={0},date={1},id={2}'.format(pre_proj,
                                                                                          pre_date,
                                                                                          pre_num))

            # Check prearchive status, must be Ready to archive
            if pre_stat != 'READY':
                LOGGER.warn('session not ready, status=' + pre_stat)
                continue

            # Look for a matching scan in RedCap
            rc_scan = None
            for r in self.rc_list:
                vuiis_label = str(r[self.rc_vuiis_field])
                rc_proj = str(r[self.rc_proj_field])
                if self.proj_map == '1':
                    rc_proj = PROJ_DICT[rc_proj]

                # Remove PI prefix if found
                if re.match("[A-Za-z]+_\d+", vuiis_label):
                    vuiis_label = vuiis_label.split('_')[1]

                if vuiis_label != pre_num:
                    continue
                elif r[self.rc_date_field] != pre_date:
                    msg = 'date mismatch:{0}, REDCap={1}, XNAT={2}'.format(pre_num,
                                                                           r[self.rc_date_field],
                                                                           pre_date)
                    self.report_error(msg)
                    continue
                elif rc_proj != pre_proj:
                    msg = 'project mismatch:{0}, REDCap={1}, XNAT={2}'.format(pre_num, rc_proj,
                                                                              pre_proj)

                    self.report_error(msg)
                    continue
                else:
                    rc_scan = r
                    break

            if not rc_scan:
                msg = 'no match found for Prearchive Session:' + pre_num
                self.report_error(msg)
                continue

            LOGGER.debug('found match for Prearchive Session:' + pre_num)
            for f in self.rc_fields:
                LOGGER.debug(f + '=' + str(rc_scan[f]))

            proj_label = str(rc_scan[self.rc_proj_field])
            if self.proj_map == '1':
                proj_label = PROJ_DICT[proj_label]

            subj_label = str(rc_scan[self.rc_subj_field])
            exp_label = str(rc_scan[self.rc_sess_field])
            vuiis_label = str(rc_scan[self.rc_vuiis_field])

            # Remove PI prefix if found
            if re.match("[A-Za-z]+_\d+", vuiis_label):
                vuiis_label = vuiis_label.split('_')[1]

            if not exp_label.startswith(subj_label):
                msg = 'invalid session ID, must use subject ID as prefix:' + exp_label
                self.report_error(msg)
                continue

            # Check for project, must exist
            xnat_proj = self.xnat.select.project(proj_label)
            if not xnat_proj.exists():
                msg = 'project does not exist in XNAT:' + proj_label
                self.report_error(msg)
                continue

            # Check for exp, cannot exist
            xnat_subj = xnat_proj.subject(subj_label)
            if xnat_subj.exists():
                xnat_exp = xnat_subj.experiment(exp_label)
                if xnat_exp.exists():
                    msg = 'session already exists in XNAT:'
                    msg += proj_label + ':' + subj_label + ':' + exp_label
                    self.report_error(msg)
                    continue

            # Archive it
            LOGGER.info('archiving Session:{0} to {1}/{2}/{3}...'.format(pre_uri,
                                                                         proj_label,
                                                                         subj_label,
                                                                         exp_label))

            result = self.archive_prearchive_session(pre_uri, proj_label, subj_label, exp_label)
            LOGGER.debug('RESULT=' + result)

    def crosscheck_redcap(self):
        """ cross check between XNAT and REDCap """
        # Check scans in REDCap against XNAT Archive
        LOGGER.info('checking RedCap against XNAT Archive...')
        for r in self.rc_list:
            exp_label = str(r[self.rc_sess_field])

            # Check for blanks in RC
            if not r[self.rc_subj_field]:
                msg = 'blank subject:session=' + exp_label
                self.report_error(msg)
                continue
            elif not r[self.rc_vuiis_field]:
                msg = 'blank vuiis number, session=' + exp_label
                self.report_error(msg)
                continue

            proj_label = str(r[self.rc_proj_field])
            if self.proj_map == '1':
                proj_label = PROJ_DICT[proj_label]

            subj_label = str(r[self.rc_subj_field])
            vuiis_label = str(r[self.rc_vuiis_field])

            # Check for project
            xnat_proj = self.xnat.select.project(proj_label)
            if not xnat_proj.exists():
                msg = 'project specified in REDCap does not exist in XNAT:'
                msg += 'vuiis#={0}, sess#={1}, subj#={2}, project={3}'.format(vuiis_label,
                                                                              exp_label,
                                                                              subj_label,
                                                                              proj_label)
                self.report_error(msg)
                continue

            # Check for subject
            xnat_subj = xnat_proj.subject(subj_label)
            if not xnat_subj.exists():
                msg = 'subject specified in REDCap does not exist in XNAT:'
                msg += 'vuiis#={0}, sess#={1}, subj#={2}, project={3}'.format(vuiis_label,
                                                                              exp_label,
                                                                              subj_label,
                                                                              proj_label)

                self.report_error(msg)
                continue

            # Check for exp
            xnat_exp = xnat_subj.experiment(exp_label)
            if not xnat_exp.exists():
                msg = 'session specified in REDCap does not exist in XNAT:'
                msg += 'vuiis#={0}, sess#={1}, project={2}'.format(vuiis_label,
                                                                   exp_label,
                                                                   proj_label)
                self.report_error(msg)
                continue

            LOGGER.debug('RedCap session:scan={0}, already archived.'.format(exp_label))

    def report_error(self, message):
        """ send message for report """
        LOGGER.error(message)
        self.report('''ERROR: {message}'''.format(message=message))

    def prerun(self, settings_filename=''):
        """ prerun function overridden from base-class """
        LOGGER.info('loading XNAT...')
        self.xnat = XnatUtils.get_interface()

        LOGGER.info('loading REDCap...')
        self.load_redcap()

        LOGGER.info('checking Projects of Prearchive Scans..')
        self.load_prearchive()
        self.check_projects()

        LOGGER.info('do archiving...')
        self.load_prearchive()
        self.do_archiving()

        LOGGER.info('crosscheck redcap...')
        self.crosscheck_redcap()

        self.xnat.disconnect()

        if self.send_an_email:
            LOGGER.info('sending email with report of errors')
            self.send_report()

        LOGGER.info('DONE')

    def afterrun(self, xnat, project):
        """ afterrun function overridden from base-class """
        pass

    def needs_run(self, csess, xnat):
        """ needs_run function to check if it needs to run on the session (never) """
        return False

    def run(self, sess_info, sess_obj):
        """ run main function called per session """
        pass
