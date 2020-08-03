from __future__ import division
import sys, os, copy, time, datetime
from multiprocessing import Process,cpu_count
from multiprocessing import Pool
import subprocess
# import glob, inspect
import argparse, pipes, math, psutil
try:
    resource
except Exception as e:
    pass

def GetMemoryUsage(print_status = False):
    available_mem_gib = None
    available_mem_percent = None
    try:
        available_mem_gib = psutil.virtual_memory().available / float (2 ** 30)
        available_mem_percent = 100 - psutil.virtual_memory().percent
    except:
        pass

    used_mem_gib = None
    # used_mem_gib2 = None
    try:
        process = psutil.Process(os.getpid())
        used_mem_gib = process.memory_info().rss / float(2 ** 30)
    except:
        try:
            used_mem_gib = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / float(2**20)
        except: pass

    if used_mem_gib != None and available_mem_gib != None:
        if print_status:
            print('memory used by the app: %g Gb'%used_mem_gib)
            print('memory available: %g Gb'%available_mem_gib)
        return used_mem_gib, available_mem_gib, available_mem_percent
    else:
        if print_status:
            print('Cannot identify sizes of total, available and used memory.')
        return None, None, None

def GetCPUUsage(print_status = False):
    try:
        CPU_usage = psutil.cpu_percent()
        if print_status:
            print('CPU usage: %.1f%% (all cores)', CPU_usage)
    except:
        CPU_usage = None

    return CPU_usage

def GetTime():
    ts = time.time()
    return  datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')

def Run_Single_CommandLine(CommandLine):
    code = os.system(CommandLine)
    return code

def Run_Single_CommandLine_And_Check_ErrorCode(CommandLine):
    code = os.system(CommandLine)
    if not (code in [0]):
        with open('Run.CL.multithread.log','a') as input:
            input.write('Command "%s" FAILED with error code %d'%(CommandLine,code) + '\n')
        # print('Command "%s" FAILED with error code %d'%(CommandLine,code))
    return code

def Run_CommandLineList_And_Check_ErrorCode(CommandLine_List):
    for CommandLine in CommandLine_List:
        code = os.system(CommandLine)
        if code != 0:
            with open('Run.CL.multithread.log','a') as input:
                input.write('Command "%s" FAILED with error code %d'%(CommandLine,code) + '\n')
            # print('Command "%s" FAILED with error code %d'%(CommandLine,code))
            return code
    return 0

def Run_Single_CommandLine___v3(arguments, stdout_output_file = None, stderr_output_file = None, append_to_logs = False):
    if stderr_output_file is None:  stderr_output_file = stdout_output_file
    
    if stdout_output_file is not None:
        if not os.path.exists(stdout_output_file):  f_stdout = open(stdout_output_file, 'w', buffering = 1)
        else:   f_stdout = open(stdout_output_file, ('a' if append_to_logs else 'w'), buffering = 1)
    else:
        f_stdout = subprocess.PIPE

    if stderr_output_file == stdout_output_file:
        f_stderr = f_stdout
    elif stderr_output_file is not None:
        if not os.path.exists(stderr_output_file):  f_stderr = open(stderr_output_file, 'w', buffering = 1)
        else:   f_stderr = open(stderr_output_file, ('a' if append_to_logs else 'w'), buffering = 1)
    else:
        f_stderr = subprocess.PIPE

    return(subprocess.run(arguments, stderr = f_stderr, stdout = f_stdout))



def Run_Single_CommandLine___v2_plus_output(arguments):
    return(subprocess.run(arguments))


def Run_Single_CommandLine_And_CheckFile(Pars):
    CommandLine,FileName = Pars
    code = os.system(CommandLine)
    if code != 0:
        with open('Run.CL.multithread.log','a') as input:
            input.write('Command "%s" FAILED with error code %d'%(CommandLine,code) + '\n')
        print('Command "%s" FAILED with error code %d'%(CommandLine,code))
        return code
    if not os.path.exists(FileName):
        print('Command "%s" FAILED'%CommandLine)
        print('Output file %s does not exist'%FileName)
        return code
    if os.path.getsize(FileName) <= 0:
        with open('Run.CL.multithread.log','a') as input:
            input.write('Command "%s" FAILED'%CommandLine + 'Output file %s has zero size'%FileName + '\n')
        print('Command "%s" FAILED'%CommandLine + 'Output file %s has zero size'%FileName)
        return code
        # exit()
    return code

def Run_CommandLines_In_MultiThread(CommandLine_List,ThreadsCount):
    pool = Pool(processes = ThreadsCount)
    pool.map(Run_Single_CommandLine, CommandLine_List)
    pool.close()
    pool.join()

def Run_CommandLines_In_MultiThread_And_Check_ErrorCode(CommandLine_List,ThreadsCount):
    pool = Pool(processes = ThreadsCount)
    Codes = pool.map(Run_Single_CommandLine_And_Check_ErrorCode, CommandLine_List)
    pool.close()
    pool.join()
    errors_count = 0
    for i in range(len(CommandLine_List)):
        if Codes[i] != 0:
            print('\n[Summary] Command "%s" FAILED with error code %d'%(str(CommandLine_List[i]),Codes[i]))
            errors_count += 1
    if errors_count > 0:
        print('\nTotal %d of %d commands FAILED\n'%(errors_count,len(CommandLine_List)))
        exit(127)

def Run_CommandLines_ListPack_In_MultiThread_And_Check_ErrorCode(CommandLine_ListPack,ThreadsCount):
    pool = Pool(processes = ThreadsCount)
    Codes = pool.map(Run_CommandLineList_And_Check_ErrorCode, CommandLine_ListPack)
    pool.close()
    pool.join()
    errors_count = 0
    for i in range(len(CommandLine_ListPack)):
        if Codes[i] != 0:
            print('[Summary] Command "%s" FAILED with error code %d'%(str(CommandLine_ListPack[i]),Codes[i]))
            errors_count += 1
    if errors_count > 0:
        print('\nTotal %d of %d commands FAILED\n'%(errors_count,len(CommandLine_ListPack)))
        exit(127)

def Run_CommandLines_In_MultiThread_And_CheckOutputFiles(CommandLine_List,OutputFiles_List,ThreadsCount):
    pool = Pool(processes = ThreadsCount)
    pool.map(Run_Single_CommandLine_And_CheckFile, list(zip(CommandLine_List,OutputFiles_List)))
    pool.close()
    pool.join()


def Log_to_File(message, LogFile_name = None, show_in_console = True, timestamp = True):
    # if LogFile_name is None:
    #     LogFile_name = './Meta-R_log.txt'

    if LogFile_name is not None:
        try:
            LogFile = open(LogFile_name, 'a')
        except:
            print('Cannot open log file ' + LogFile_name)
            return

    if timestamp:
        message = '[%s] '%GetTime() + str(message)
    
    if show_in_console:
        print(str(message))

    if LogFile_name is not None:
        LogFile.write(str(message) + '\n')
        LogFile.close()

def Run_CommandLines_In_MultiThread___v3(arguments_List, ThreadsCount, main_LogFile_name = None,
                                         individual_LogFile_names = None, individual_err_LogFile_names = 'same', append_logs = False, 
                                         suppress_logging_to_console = False, command_names = None,
                                         terminal_columns = None, touch_ok_file = False, ok_files = []):

    pool = Pool(processes = ThreadsCount)
    # rs = pool.imap_unordered(Run_Single_CommandLine___v2, arguments_List)

    command_names_present = True
    command_lines = [' '.join(['"%s"'%x for x in arguments]) for arguments in arguments_List]
    if command_names is None:
        command_names = [' '.join(['"%s"'%x for x in arguments]) for arguments in arguments_List]
        command_names_present = False
    elif len(command_names) != len(arguments_List):
        raise('length of arguments_List (%d) should be equal to the length of command_names (%d)'%(len(arguments_List), len(command_names)))

    if individual_LogFile_names is None:
        individual_LogFile_names = [None]*len(arguments_List)
        individual_LogFile_names_present = False
    elif len(individual_LogFile_names) != len(arguments_List):
        raise('length of arguments_List (%d) should be equal to the length of individual_LogFile_names (%d)'%(len(arguments_List), len(individual_LogFile_names)))
    else:
        individual_LogFile_names_present = True

    if individual_err_LogFile_names is not None:
        if individual_err_LogFile_names == 'same':
            individual_err_LogFile_names = individual_LogFile_names
            individual_err_LogFile_names_present = individual_LogFile_names_present
        elif len(individual_err_LogFile_names) != len(arguments_List):
            raise('length of arguments_List (%d) should be equal to the length of individual_err_LogFile_names (%d)'%(len(arguments_List), len(individual_err_LogFile_names)))
        else:
            individual_err_LogFile_names_present = True
    else:
        individual_err_LogFile_names = [None]*len(arguments_List)
        individual_err_LogFile_names_present = False

    stdout_and_stderr_log_files_are_same = (individual_err_LogFile_names_present and individual_LogFile_names_present and str(individual_err_LogFile_names) == str(individual_LogFile_names))

    append_logs = [append_logs]*len(arguments_List)

    if touch_ok_file:
        if len(ok_files) != len(arguments_List):
            raise('length of arguments_List (%d) should be equal to the length of ok_files (%d)'%(len(arguments_List), len(ok_files)))


    pool = Pool(processes = ThreadsCount)
    all_results = pool.starmap_async(Run_Single_CommandLine___v3, zip(arguments_List, individual_LogFile_names, individual_err_LogFile_names, append_logs))
    pool.close() # No more work
    print('\n')
    prev_text = 'None'
    while (True):
        if (all_results.ready()): break
        remaining = all_results._number_left
        # print("Waiting for " + str(remaining) + " tasks to complete...")
        used_mem_gib, available_mem_gib, available_mem_percent = GetMemoryUsage()
        cpu_percent = GetCPUUsage()
        if available_mem_gib is not None and cpu_percent is not None:
            msg = '%d of %d tasks completed. %.1f Gb RAM available (%.0f%%). CPU load (all cores): %.1f%%'%((len(arguments_List) - remaining), len(arguments_List), available_mem_gib, available_mem_percent, cpu_percent)

            # msg = '%d of %d tasks completed. %s Gb RAM available (%s), %s / %s Gb used'%((len(arguments_List) - remaining), len(arguments_List), str(available_mem_gib), str(available_mem_percent), str(used_mem_gib), str(used_mem_gib2))
        else:
            msg = '%d of %d tasks completed'%((len(arguments_List) - remaining), len(arguments_List))

        if msg != prev_text:
            sys.stdout.write('\r   ' + msg + '           \r')
            Log_to_File(msg, main_LogFile_name, show_in_console = False)
            prev_text = msg
        # print(all_results.tasks_remaining())
        time.sleep(0.5)

    # pool.close() # No more work
    pool.join() # Wait for completion
    if terminal_columns is None:
        terminal_rows = 40
        terminal_columns = 120
        try:
            terminal_rows, terminal_columns = os.popen('stty size', 'r').read().split()
            terminal_columns = int(terminal_columns)
            terminal_rows = int(terminal_rows)
        except:  pass

    failed_results = []
    failed_results_n = []
    succeed_results = []
    succeed_results_n = []
    all_results__get = all_results.get(timeout=1)
    # for r in 
    # print(all_results__get)
    for nr in range(len(all_results__get)):
        result = all_results__get[nr]
        if result.returncode == 0:
            succeed_results += [result]
            succeed_results_n += [nr]
            if touch_ok_file:
                os.system('touch "%s"'%ok_files[nr])
        else:
            failed_results += [result]
            failed_results_n += [nr]

    if len(failed_results) == 0:
        msg = 'All %d tasks completed successfully'%len(arguments_List)
        Log_to_File('\r%s'%msg + ' '*max(0, terminal_columns - len(msg)), main_LogFile_name)
    else:
        Log_to_File('\r' + '*'*terminal_columns, main_LogFile_name, timestamp = False)
        Log_to_File('Warning. %d of %d tasks failed:\n%s'%(len(failed_results), len(arguments_List),
            '\n'.join([' '*8 + str(command_names[x]) for x in failed_results_n])), main_LogFile_name)

    
    if suppress_logging_to_console == 'if_fails' and len(failed_results) > 0:
        show_in_console_failed = True
        show_in_console_success = False
    else:
        show_in_console_failed = not suppress_logging_to_console
        show_in_console_success = not suppress_logging_to_console

    Log_to_File('\n\n Successful commands logs:\n\n', main_LogFile_name, show_in_console = show_in_console_success, timestamp = False)
    for nr in succeed_results_n:
        Log_to_File("%d. %s"%(nr + 1, command_names[nr]), main_LogFile_name, show_in_console = show_in_console_success, timestamp = False)
        if command_names_present:
            Log_to_File("Command line: " + command_lines[nr], main_LogFile_name, show_in_console = show_in_console_success, timestamp = False)

        if stdout_and_stderr_log_files_are_same:
            Log_to_File("\nOutput:", main_LogFile_name, show_in_console = show_in_console_success, timestamp = False)
            Log_to_File(open(individual_LogFile_names[nr], 'r').read(), main_LogFile_name, show_in_console = show_in_console_success, timestamp = False)
        else:
            if individual_err_LogFile_names_present:
                stderr_record = open(individual_err_LogFile_names[nr], 'r').read()
            else:
                stderr_record = all_results__get[nr].stderr.decode('utf-8')
            Log_to_File("\nSTDERR:", main_LogFile_name, show_in_console = show_in_console_success, timestamp = False)
            Log_to_File(stderr_record, main_LogFile_name, show_in_console = show_in_console_success, timestamp = False)

            if individual_LogFile_names_present:
                stdout_record = open(individual_LogFile_names[nr], 'r').read()
            else:
                stdout_record = all_results__get[nr].stdout.decode('utf-8')
            Log_to_File("STDOUT:", main_LogFile_name, show_in_console = show_in_console_success, timestamp = False)
            Log_to_File(stdout_record + '\n', main_LogFile_name, show_in_console = show_in_console_success, timestamp = False)

        Log_to_File('*'*terminal_columns + '\n', show_in_console = show_in_console_success, timestamp = False)

    if len(failed_results) > 0:
        Log_to_File('\n\n Failed commands logs:\n\n', main_LogFile_name, show_in_console = show_in_console_failed)
        for nr in failed_results_n:
            Log_to_File("%d. %s"%(nr + 1, command_names[nr]), main_LogFile_name, show_in_console = show_in_console_failed, timestamp = False)
            if command_names_present:
                Log_to_File("Command line: " + command_lines[nr], main_LogFile_name, show_in_console = show_in_console_failed, timestamp = False)
            
            if stdout_and_stderr_log_files_are_same:
                Log_to_File("\nOutput:", main_LogFile_name, show_in_console = show_in_console_failed, timestamp = False)
                Log_to_File(open(individual_LogFile_names[nr], 'r').read(), main_LogFile_name, show_in_console = show_in_console_failed, timestamp = False)
            else:
                if individual_err_LogFile_names_present:
                    stderr_record = open(individual_err_LogFile_names[nr], 'r').read()
                else:
                    stderr_record = all_results__get[nr].stderr.decode('utf-8')
                Log_to_File("\nSTDERR:", main_LogFile_name, show_in_console = show_in_console_failed, timestamp = False)
                Log_to_File(stderr_record, main_LogFile_name, show_in_console = show_in_console_failed, timestamp = False)
                
                if individual_LogFile_names_present:
                    stdout_record = open(individual_LogFile_names[nr], 'r').read()
                else:
                    stdout_record = all_results__get[nr].stdout.decode('utf-8')
                Log_to_File("STDOUT:", main_LogFile_name, show_in_console = show_in_console_failed, timestamp = False)
                Log_to_File(stdout_record + '\n', main_LogFile_name, show_in_console = show_in_console_failed, timestamp = False)

            Log_to_File('*'*terminal_columns + '\n', show_in_console = show_in_console_success, timestamp = False)

        Log_to_File('Finally, %d of %d tasks failed:\n%s\n'%(len(failed_results), len(arguments_List),
            '\n'.join([' '*7 + '%d. '%(xn + 1) + str(command_names[failed_results_n[xn]]) for xn in range(len(failed_results_n))])), main_LogFile_name, show_in_console = show_in_console_failed)
    else:
        Log_to_File('\rFinally all %d tasks completed successfully                                \n'%(len(arguments_List)), main_LogFile_name, show_in_console = show_in_console_success)
                                                                   

    return(all_results)

    # for i, _ in enumerate(p.imap_unordered(Run_Single_CommandLine___v2, arguments_List), 1):
    #     sys.stderr.write('\rdone {0:%}'.format(i/len(arguments_List)))

    # pool = Pool(processes = ThreadsCount)
    # pool.map(Run_Single_CommandLine___v2, arguments_List)
    # pool.close()
    # pool.join()

def ToBool(string):
    string = string.casefold()
    if string in ['yes','y','on','enable','true','da','t']:
        return True
    elif string in ['no','n','off','disable','false','net','f']:
        return False
    print('Parameter "%s" should be either "yes" or "no". Please correct'%string)
    exit()



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-cf', '--cl-file', dest = 'command_list_filename', nargs='?', action='store', required = False, default = None, help = '')
    parser.add_argument('-cd', '--commands', dest = 'commands', nargs='*', action='store', required = False, default = None, help = '')
    parser.add_argument('-logs', '--log-files', dest = 'log_files', nargs='*', action='store', required = False, default = None, help = '')
    parser.add_argument('-err-logs', '--err-log-files', dest = 'err_log_files', nargs='*', action='store', required = False, default = 'same_as_log_files', help = '')
    parser.add_argument('-main-log', '--main-log-file', dest = 'main_log_file', nargs='?', action='store', required = False, default = None, help = '')
    parser.add_argument('-a', '--append-logs', dest = 'append_logs', nargs='?', action='store', required = False, default = 'no', help = '')
    parser.add_argument('-sup', '--suppress-logging-to-console', dest = 'suppress_logging_to_console', nargs='?', action='store', required = False, default = 'no', help = '')
    parser.add_argument('--sep', dest = 'arguments_separator', nargs='?', action='store', required=False, default = ' ', help = '')
    parser.add_argument('-t', '--threads', dest = 'threads', nargs='?', action='store', required=False, default = None, help = '')
    if '-h' in sys.argv or '--help' in sys.argv or len(sys.argv) == 1:
        print('Usage:\n')
        print('python3 Run.CL.multithread.py -cd "blastx aaa bbb" "blastp ccc ddd" -t 16')
        exit(0)
    args = parser.parse_args()

    args.append_logs = ToBool(args.append_logs)
    if args.suppress_logging_to_console != 'if_fails':
        args.suppress_logging_to_console = ToBool(args.suppress_logging_to_console)

    if (args.commands is None) == (args.command_list_filename is None):
        print('Please specify either commands or command list file (-cf or -cd arguments)')
        exit(127)

    if args.log_files is not None and args.commands is None:
        print('Logging to files is supported only for -cd argument')
        exit(127)

    if args.log_files is not None:
        if len(args.log_files) != len(args.commands):
            print('length of args.log_files (-logs) should be equal to the length of args.commands (-cd)')
            exit(127)
        
        if args.main_log_file is not None:
            print('please specify either main log file name or individual log files')
            exit(127)

    if args.err_log_files == 'same_as_log_files':  args.err_log_files = args.log_files
    elif args.err_log_files is not None and args.log_files is not None:
        if len(args.err_log_files) != len(args.log_files):
            print('The number of -logs and -err-logs arguments should be equal')
            exit(127)

    if args.threads == None: threads = cpu_count()
    else: args.threads = int(args.threads)

    if args.command_list_filename is not None:
        Run_CommandLines_In_MultiThread_And_Check_ErrorCode([x.rstrip() for x in open(args.command_list_filename).readlines()], args.threads)
    else:
        arguments_List = [[ar.strip('\'"\r\n\t ').rstrip('\'"\r\n\t ') for ar in com.split(args.arguments_separator)] for com in args.commands]
        
        Run_CommandLines_In_MultiThread___v3(arguments_List, args.threads, main_LogFile_name = args.main_log_file,
                                         individual_LogFile_names = args.log_files, individual_err_LogFile_names = args.err_log_files,
                                         append_logs = args.append_logs, suppress_logging_to_console = args.suppress_logging_to_console,
                                         command_names = None,
                                         terminal_columns = None, touch_ok_file = False, ok_files = [])

