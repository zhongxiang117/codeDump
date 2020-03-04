import os
import glob
import argparse

def file_size_check(path,fsize=20):
    """This function is used to check the file existence and size,
       the unit of size is in megabety"""

    fndx = True
    try:
        sizetmp = os.stat(path).st_size
        if sizetmp/1024/1024 > fsize:
            print('Error: the file size is far larger than %f MB' % fsize)
            print('Error path: ',path)
            fndx = False         
    except IOError:
        print('Error : cannot open the file!')
        print('Error path: ',path)
        fndx = False
    return fndx


def file_gen_new(fname,fextend='txt',foriginal=True,bool_dot=True):
    """This function is used to make new file but without overriding the
       old one. By default, the "string" after the dot-notation in the
       input file name is used as the final file extension, which means
       it will override the parameter in 'fextend', this behavior can be
       turned off by set 'bool_dot' to False. The parameter 'foriginal'
       is used to set filename be counted or not"""
    
    filename = fname
    pos = filename.rfind('.')
    if bool_dot and pos != -1:
        fname = filename[:pos]
        fextend = filename[pos:]
    else:
        fextend = '.' + fextend
    
    if foriginal is True:
        try:
            f = open(fname + fextend)
            f.close()
        except:
            return fname + fextend
        
    i = 1
    filename = fname
    while True:
        fname = filename + '_' + str(i) + fextend
        try:
            f = open(fname)
            f.close()
        except:
            break
        i += 1
    return fname


def file_decoding(file):
    """Decode the input file, and return the text file list,
       if the file cannot be successfully processed,
       a empty list is returned"""

    encodings = ('utf_8', 'utf_16', 'utf_16_le', 'utf_16_be','utf_32',
                 'utf_32_be', 'utf_32_le','cp850' , 'cp437', 'cp852',
                 'cp1252', 'cp1250' , 'ascii','utf_8_sig', 'big5',
                 'big5hkscs', 'cp037', 'cp424', 'cp500','cp720', 'cp737',
                 'cp775', 'cp855', 'cp856', 'cp857','cp858', 'cp860',
                 'cp861', 'cp862', 'cp863', 'cp864','cp865', 'cp866',
                 'cp869', 'cp874', 'cp875', 'cp932','cp949', 'cp950',
                 'cp1006', 'cp1026', 'cp1140', 'cp1251','cp1253',
                 'cp1254', 'cp1255', 'cp1256', 'cp1257', 'cp1258',
                 'euc_jp', 'euc_jis_2004', 'euc_jisx0213','euc_kr',
                 'gb2312', 'gbk', 'gb18030', 'hz', 'iso2022_jp',
                 'iso2022_jp_1', 'iso2022_jp_2', 'iso2022_jp_2004',
                 'iso2022_jp_3', 'iso2022_jp_ext', 'iso2022_kr', 'latin_1',
                 'iso8859_2', 'iso8859_3', 'iso8859_4', 'iso8859_5',
                 'iso8859_6', 'iso8859_7', 'iso8859_8', 'iso8859_9',
                 'iso8859_10', 'iso8859_13', 'iso8859_14', 'iso8859_15',
                 'iso8859_16', 'johab', 'koi8_r', 'koi8_u', 'mac_cyrillic',
                 'mac_greek', 'mac_iceland', 'mac_latin2', 'mac_roman',
                 'mac_turkish', 'ptcp154', 'shift_jis', 'shift_jis_2004',
                 'shift_jisx0213')

    f = open(file,mode='rb').read()

    profile = []
    for e in encodings:
        try:
            fdec = f.decode(e)
            if 'ame' in fdec:
                profile = fdec.split('\n')
                break
        except:
            pass

    if len(profile) == 0:
        return []
    return profile


class Average_method(object):
    def __init__(self):
        pass
        
    def average_all(self,grade):
        half = len(grade) // 2
        aver_quiz = 1.0 * sum(grade[:half]) / half
        aver_report = 1.0 * sum(grade[half:]) / half
        return aver_quiz,aver_report
                                   

    def average_drop_min_quiz_only(self,grade):
        half = len(grade) // 2
        subsum = sum(grade[:half]) - min(grade[:half])
        aver_quiz = 1.0 * subsum / (half-1)
        
        aver_report = 1.0 * sum(grade[half:]) / half
        return aver_quiz,aver_report


    def average_drop_min_report_only(self,grade):
        half = len(grade) // 2
        aver_quiz = 1.0 * sum(grade[:half]) / half
        
        subsum = sum(grade[half:]) - min(grade[half:])
        aver_report = 1.0 * subsum / (half-1)
        return aver_quiz,aver_report


    def average_drop_min_experiment_basedOnQuiz(self,grade):
        half = len(grade) // 2
        gmin = min(grade[:half])
        ndx = grade[:half].index(gmin)

        subsum = sum(grade[:half]) - gmin
        aver_quiz = 1.0 * subsum / (half-1)
        
        gdrop = grade[ndx+half]
        subsum = sum(grade[half:]) - gdrop
        aver_report = 1.0 * subsum / (half-1)
        return aver_quiz,aver_report


    def average_drop_min_experiment_basedOnReport(self,grade):
        half = len(grade) // 2
        gmin = min(grade[half:])
        ndx = grade[half:].index(gmin)

        subsum = sum(grade[half:]) - gmin
        aver_report = 1.0 * subsum / (half-1)
        
        gdrop = grade[ndx]
        subsum = sum(grade[:half]) - gdrop
        aver_quiz = 1.0 * subsum / (half-1)
        return aver_quiz,aver_report


    def average_drop_min_experiment(self,grade,ratio=[0.20,0.55]):
        half = len(grade) // 2
        reflist = []
        for i in range(half):
            reflist.append(grade[:half][i]*ratio[0] + grade[half:][i]*ratio[1])
        ndx = reflist.index(min(reflist))

        subsum = sum(grade[:half]) - grade[ndx]
        aver_quiz = 1.0 * subsum / (half-1)

        subsum = sum(grade[half:]) - grade[half:][ndx]
        aver_report = 1.0 * subsum / (half-1)
        return aver_quiz,aver_report        


def file_process(file):   
    """
    Return Value:
        gradelist, studentlist

    
    Structure:
        gradelist:
            QuizGrade..., ReportGrade..., EOSquiz
        studentlist:
            lastName, firstName, C_Number,
    """
    
    profile = []
    experiment = []
    style_ratio = []
    for line in file:
        ltmp = line.split('\t')
        if len(ltmp) > 3:
            if ltmp[0].lower().find('name') != -1 and \
               ltmp[0].lower().find('last') != -1:
                for ndx in range(len(ltmp)):
                    if ltmp[ndx].lower().find('q:') != -1 or \
                       ltmp[ndx].lower().find('r:') != -1 or \
                       ltmp[ndx].lower().find('end') != -1:
                        experiment.append(ndx)

                if len(experiment) == 0 or len(experiment) % 2 != 1: 
                    print('Error: the name line entry is wrong')
                    print(line)
                    exit()

                line = line.lower()
                if line.find('style') != -1:
                    st = line[line.find('style'):].split(',')
                    for tmp in st:
                        if tmp.find('style') != -1:
                            try:
                                st = tmp.split(':')[1]
                                st = st.split('\'')[0].split('"')[0]
                                nm = int(st)                                
                            except:
                                nm = 0
                            style_ratio.insert(0,nm)
                            break
                if line.find('ratio') != -1:
                    st = line[line.find('ratio'):].split(',')
                    for tmp in st:
                        if tmp.find('ratio') != -1:
                            riostr = tmp.split(':')
                            if len(riostr) >= 2:
                                rio = riostr[1].split('-')
                                if len(rio) == 3:
                                    try:
                                        rio_q = float(rio[0].split()[0])
                                        rio_r = float(rio[1].split()[0])

                                        f = rio[2].split()[0]
                                        f = f.split('\'')[0].split('"')[0]
                                        rio_f = float(f)

                                        style_ratio.append([rio_q,rio_r,rio_f])
                                    except:
                                        style_ratio.append([])                        
            elif len(experiment) != 0 and len(ltmp) > max(experiment):
                profile.append(ltmp)

    gradelist = []
    studentlist = []
    for student in profile:
        ls = []
        for ndx in experiment:
            try:
                tmpstr = ''
                for ch in student[ndx]:
                    if ch != '\'' and ch != '"' and ch != '\n':
                        tmpstr += ch
                t = int(tmpstr.split()[0])
            except ValueError:
                t = 0
            ls.append(t)
        gradelist.append(ls)
        studentlist.append([student[0],student[1],student[2]])
    return gradelist,studentlist,style_ratio



parser = argparse.ArgumentParser(description='Chemistry Lab Grades Process',allow_abbrev=False)
parser.add_argument('-f','--gradefile',help='Lab grades file path')
parser.add_argument('-r','--ratio',help='Ratio amonge Quiz, Reports, and EndOfSemesterQuiz,\
                     only three arguments are accepted, default is 0.20, 0.55, 0.25',nargs=3,type=float)
parser.add_argument('-t','--style',help='Printout style, 0 prints all, 1 prints "Average all the \
                     grades", 2 prints "Average by dropping one lowest quiz", 3 prints "Average \
                     by dropping one lowest report", 4 prints "Average by dropping one experiment \
                     based on quiz", 5 prints "Average by dropping one experiment based on report", \
                     6 prints "Average by dropping one lowest experiment based on course grade \
                     ratio/scheme", 7 prints "Average of the best performance". Default is 0, \
                     prints all',nargs=1,type=int)
parser.add_argument('-o','--fname',help='Output file name, default is LabProcessResults')

args = parser.parse_args()

# Note, the sequence is Quiz, Report, EOSquiz
if args.ratio is not None:
    if sum(args.ratio) - 1.0 > 0.0001:
        print('Warning: the ratio is not based on 100 percent!')
        print('Do you want to continue? y/yes, else quit',end='    ')
        getstr = input()
        if getstr.lower() != 'y' and getstr.lower() != 'yes':
            print('Warning: you have decided to quit!')
            print('Exiting ...')
            exit()
    ratio = args.ratio

    
if args.gradefile is None:
    flist = glob.glob('CHM*')
else:
    flist = [args.gradefile]

# Instance Average Method
AM = Average_method()


def gradeCalc(aver_quiz,aver_report,eos,ratio):
    final = aver_quiz*ratio[0]*100/25 + aver_report*ratio[1]*100/50 + eos*ratio[2]
    line = '{:>7.2f}\t{:>7.2f}\t{:>5.0f}\t{:>7.2f}'.format(aver_quiz,aver_report,eos,final)
    return final, line

outlist = []
stylelist = []
title = 'Last Name         \tFirst Name        \t     C Number\t   Quiz\t Report\t  EOS\t  Final\n'
for file in flist:
    score = []
    outfile = []
    if file_size_check(file):

        profile = file_decoding(file)
        if len(profile) == 0:
            continue
        else:
            gradelist,std,style_ratio = file_process(profile)
            if len(std) == 0:
                continue
            
        if args.ratio is None:
            if len(style_ratio) == 0:
                ratio = [0.20,0.55,0.25]
            elif len(style_ratio) == 1:
                if isinstance(style_ratio[0],int):
                    ratio = [0.20,0.55,0.25]
                else:
                    if len(style_ratio[0]) == 0:
                        ratio = [0.20,0.55,0.25]
                    else:
                        ratio = style_ratio[0]
            else:
                if isinstance(style_ratio[0],list):
                    tmp = style_ratio[0]
                else:
                    tmp = style_ratio[1]
                    
                ratio = [0.20,0.55,0.25] if len(tmp) == 0 else tmp

        if args.style is None:
            if len(style_ratio) == 0:
                style = 0
            elif len(style_ratio) == 1:
                if isinstance(style_ratio,int):
                    if style_ratio[0] in [0,1,2,3,4,5,6,7]:
                        style = style_ratio[0]
                    else:
                        style = 0
                else:
                    style = 0
            else:
                if style_ratio[0] in [0,1,2,3,4,5,6,7]:
                    style = style_ratio[0]
                else:
                    style = 0
        else:
            if args.style in [0,1,2,3,4,5,6,7]:
                style = args.style
            else:
                style = 0
        stylelist.append(style)
                    
        lp = []
        lp.append('For grade file : < %s >, the process results are;\n\n\n' % file)
        lp.append('Note: the experiment grade ratio/scheme used is:\n\n')
        lp.append('    Quiz:                  {:}\n'.format(ratio[0]))
        lp.append('    Report:                {:}\n'.format(ratio[1]))
        lp.append('    End of semester quiz:  {:}\n\n\n'.format(ratio[2]))
        outfile.append(lp)

        lp = []
        ls = []
        lp.append('Average all the grades:\n\n')
        lp.append(title)
        for nm in range(len(gradelist)):
            aver_quiz,aver_report = AM.average_all(gradelist[nm][:-1])
            
            #line = '{:18}\t{:18}\t{:>13}\t{:>7.2f}\t{:>7.2f}\t{:>5.0f}\t{:>7.2f}'.format(std[nm][0],\
            #                   std[nm][1],std[nm][2],aver_quiz,aver_report,gradelist[nm][-1],final)

            stdinfo = '{:18}\t{:18}\t{:>13}'.format(std[nm][0],std[nm][1],std[nm][2])

            final, line = gradeCalc(aver_quiz,aver_report,gradelist[nm][-1],ratio)

            outline = stdinfo + '\t' + line + '\n'
            ls.append(final)
            lp.append(outline)
        lp.append('\n\n')
        score.append(ls)
        outfile.append(lp)


        lp = []
        ls =[]
        lp.append('Average by dropping one lowest quiz grade:\n\n')
        lp.append(title)
        for nm in range(len(gradelist)):
            aver_quiz,aver_report = AM.average_drop_min_quiz_only(gradelist[nm][:-1])
            
            stdinfo = '{:18}\t{:18}\t{:>13}'.format(std[nm][0],std[nm][1],std[nm][2])
            final, line = gradeCalc(aver_quiz,aver_report,gradelist[nm][-1],ratio)
            outline = stdinfo + '\t' + line + '\n'
            ls.append(final)
            lp.append(outline)
        lp.append('\n\n')
        score.append(ls)
        outfile.append(lp)
            

        lp = []
        ls = []
        lp.append('Average by dropping one lowest report grade:\n\n')
        lp.append(title)
        for nm in range(len(gradelist)):
            aver_quiz,aver_report = AM.average_drop_min_report_only(gradelist[nm][:-1])
            
            stdinfo = '{:18}\t{:18}\t{:>13}'.format(std[nm][0],std[nm][1],std[nm][2])
            final, line = gradeCalc(aver_quiz,aver_report,gradelist[nm][-1],ratio)
            outline = stdinfo + '\t' + line + '\n'
            ls.append(final)
            lp.append(outline)
        lp.append('\n\n')
        score.append(ls)
        outfile.append(lp)
            

        lp = []
        ls = []
        lp.append('Average by dropping one experiment based on lowest quiz:\n\n')
        lp.append(title)
        for nm in range(len(gradelist)):
            aver_quiz,aver_report = AM.average_drop_min_experiment_basedOnQuiz(gradelist[nm][:-1])
            
            stdinfo = '{:18}\t{:18}\t{:>13}'.format(std[nm][0],std[nm][1],std[nm][2])
            final, line = gradeCalc(aver_quiz,aver_report,gradelist[nm][-1],ratio)
            outline = stdinfo + '\t' + line + '\n'
            ls.append(final)
            lp.append(outline)
        lp.append('\n\n')
        score.append(ls)
        outfile.append(lp)


        lp = []
        ls = []
        lp.append('Average by dropping one experiment based on lowest report:\n\n')
        lp.append(title)
        for nm in range(len(gradelist)):
            aver_quiz,aver_report = AM.average_drop_min_experiment_basedOnReport(gradelist[nm][:-1])
            
            stdinfo = '{:18}\t{:18}\t{:>13}'.format(std[nm][0],std[nm][1],std[nm][2])
            final, line = gradeCalc(aver_quiz,aver_report,gradelist[nm][-1],ratio)
            outline = stdinfo + '\t' + line + '\n'
            ls.append(final)
            lp.append(outline)
        lp.append('\n\n')
        score.append(ls)
        outfile.append(lp)


        lp = []
        ls = []
        lp.append('Average by dropping one lowest experiment based on course grade ratio/scheme:\n\n')
        lp.append(title)
        for nm in range(len(gradelist)):
            aver_quiz,aver_report = AM.average_drop_min_experiment(gradelist[nm][:-1],ratio[:2])
            
            stdinfo = '{:18}\t{:18}\t{:>13}'.format(std[nm][0],std[nm][1],std[nm][2])
            final, line = gradeCalc(aver_quiz,aver_report,gradelist[nm][-1],ratio)
            outline = stdinfo + '\t' + line + '\n'
            ls.append(final)
            lp.append(outline)
        lp.append('\n\n')
        score.append(ls)
        outfile.append(lp)


        bestscore = []
        for ndx in range(len(score[0])):
            fgrade = []
            for stdlist in score:
                fgrade.append(stdlist[ndx])
            bestscore.append(max(fgrade))

        lp = []
        lp.append('Average of the best performance:\n\n')
        lp.append('Last Name         \tFirst Name        \t     C Number\t  Best Final\n')
        for nm in range(len(gradelist)):
            outline = '{:18}\t{:18}\t{:>13}\t{:>10.2f}\n'.format(std[nm][0],std[nm][1],std[nm][2],bestscore[nm])
            lp.append(outline)
        outfile.append(lp)
            
    outlist.append(outfile)
   
    
ndx = 0
for result in outlist:

    style = stylelist[ndx]
    if args.fname is None:
        fname = '#' + flist[ndx][:flist[ndx].rfind('.')] + '_processResult.txt'
    else:
        fname = args.fname
    ndx += 1

    fname = file_gen_new(fname)
    with open(fname,mode='wt') as f:
        
        for line in result[0]:
            f.write(line)

        if style == 0:
            for tmp in result[1:]:
                for line in tmp:
                    f.write(line)
        elif style == 1:
            for line in result[1]:
                f.write(line)
        elif style == 2:
            for line in result[2]:
                f.write(line)
        elif style == 3:
            for line in result[3]:
                f.write(line)
        elif style == 4:
            for line in result[4]:
                f.write(line)
        elif style == 5:
            for line in result[5]:
                f.write(line)
        elif style == 6:
            for line in result[6]:
                f.write(line)
        else:
            for line in result[7]:
                f.write(line)
        f.write('\n\n\n')




