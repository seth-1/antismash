# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import logging
import os

import warnings
# Don't display the SearchIO experimental warning, we know this.
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from Bio import SearchIO
from subprocess import PIPE, Popen
from helperlibs.wrappers.io import TemporaryDirectory

from antismash.config.args import Config

class RunResult():
    def __init__(self, command, stdout, stderr, return_code, piped_out, piped_err):
        self.command = command
        self.stdout = str(stdout)
        self.stderr = str(stderr)
        self.return_code = return_code
        self.stdout_piped = piped_out
        self.stderr_piped = piped_err

    def __getattr__(self, attr):
        if attr == 'stdout' and not self.stdout_piped:
            raise ValueError("stdout was redirected to file, unable to access")
        if attr == 'stderr' and not self.stderr_piped:
            raise ValueError("stdout was redirected to file, unable to access")
        return self.__dict__[attr]

    def successful(self):
        return not self.return_code

    def get_command_string(self):
        return " ".join(self.command)

def execute(commands, stdin=None, stdout=PIPE, stderr=PIPE):
    "Execute commands in a system-independent manner"

    if stdin is not None:
        stdin_redir = PIPE
    else:
        stdin_redir = None

    try:
        proc = Popen(commands, stdin=stdin_redir, stdout=stdout, stderr=stdout)
        out, err = proc.communicate(input=stdin)
        return RunResult(commands, out, err, proc.returncode,
                         stdout == PIPE, stderr == PIPE)
    except OSError as err:
        logging.debug("%r %r returned %r", commands,
                      stdin[:40] if stdin is not None else None, err)
        raise

def run_hmmsearch(query_hmmfile, target_sequence, use_tempfile=False):
    "Run hmmsearch"
    config = Config()
    command = ["hmmsearch", "--cpu", str(config.cpus),
               "-o", os.devnull, # throw away the verbose output
               "--domtblout", "result.domtab",
               query_hmmfile]

    # Allow for disabling multithreading for HMMer3 calls in the command line
    if config.get('hmmer3') and 'multithreading' in config.hmmer3 and \
            not config.hmmer3.multithreading:
        command = command[0:1] + command[3:]

    with TemporaryDirectory(change=True):
        try:
            if use_tempfile:
                with open("input.fa", 'w') as handle:
                    handle.write(target_sequence)
                command.append("input.fa")
                run_result = execute(command)
            else:
                command.append('-')
                run_result = execute(command, stdin=target_sequence)
        except OSError:
            return []
        if not run_result.successful():
            logging.error('hmmsearch returned %d: %s while searching %s',
                          run_result.return_code, run_result.stderr, query_hmmfile)
            return []
        results = list(SearchIO.parse("result.domtab", 'hmmsearch3-domtab'))
        return results

def run_hmmpress(hmmfile):
    "Run hmmpress"
    command = ['hmmpress', hmmfile]
    try:
        run_result = execute(command)
        out = run_result.stdout
        err = run_result.stderr
        retcode = run_result.return_code
    except OSError as excep:
        retcode = 1
        err = str(excep)
        out = None
    return out, err, retcode
