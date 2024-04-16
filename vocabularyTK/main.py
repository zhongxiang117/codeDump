import os
import random
import tkinter as tk
from tkinter import font
from tkinter import messagebox


# version 0.10: December 25, 2020
# version 0.20: January 3, 2021; GUI, big change


def func_read_new_vocabulary(key='vocabulary'):
    """read words from newly add vocabulary file
    
    default key: vocabulary, word, words
    file format:    key-number.txt

    return:
        str, newstr :   "word  word  word  ..."
        str, hardstr:   "word  word  word  ..."
        list, fnamelist:  [ 'filename', 'filename', ...]
    """
    def readfile(file):
        new = ''
        hard = ''
        with open(file,'rt') as f:
                while True:
                    line = f.readline()
                    if len(line) == 0:
                        break
                    line = line.strip()
                    if len(line) == 0 or line[0] == '#':
                        continue
                    if len(line.split()) == 1:
                        new += line + '  '
                    else:
                        hard += line + '  '
        return new, hard

    def get_fnamelist(key):
        fnamelist = []
        for i in range(1,1000):
            file = key + '-' + str(i) + '.txt'
            if os.path.isfile(file):
                fnamelist.append(file)
        return fnamelist


    fnamelist = get_fnamelist('vocabulary') + get_fnamelist('word')
    fnamelist += get_fnamelist('words')

    if os.path.isfile(key): fnamelist.append(key)
    if os.path.isfile(key+'.txt'): fnamelist.append(key+'.txt')
    if os.path.isfile('words.txt'): fnamelist.append('words.txt')
    if os.path.isfile('words'): fnamelist.append('words')

    newstr = ''
    hardstr = ''
    for file in fnamelist:
        new, hard = readfile(file)
        newstr += new + '  '
        hardstr += hard + '  '

    return newstr, hardstr, fnamelist


def func_read_data_file():
    """read words from data file

    file name:  data.txt

    format:
        [new]
        word  word  word ...

        [good]
        word  word  word ...

        [hard]
        word  word  word ...

    return:
        str, newstr:    "word  word  word  ..."
        str, goodstr:   "word  word  word  ..."
        str, hardstr:   "word  word  word  ..."
    """
    newstr = ''
    goodstr = ''
    hardstr = ''

    file = 'data.txt'
    if os.path.isfile(file):
        with open(file,'rt') as f:
            profile = f.readlines()
        
        i = 0
        while i < len(profile):
            line = profile[i].strip()
            if len(line) == 0 or line[0] == '#':
                i += 1
                continue
            
            if line[:5] == '[new]':
                j = i + 1
                while j < len(profile):
                    line = profile[j].strip()
                    if len(line) == 0 or line[0] == '#':
                        j += 1
                        continue
                    if line[0] == '[':
                        break
                    newstr += line + '  '
                    j += 1
            elif line[:6] == '[good]':
                j = i + 1
                while j < len(profile):
                    line = profile[j].strip()
                    if len(line) == 0 or line[0] == '#':
                        j += 1
                        continue
                    if line[0] == '[':
                        break
                    goodstr += line + '  '
                    j += 1
            elif line[:6] == '[hard]':
                j = i + 1
                while j < len(profile):
                    line = profile[j].strip()
                    if len(line) == 0 or line[0] == '#':
                        j += 1
                        continue
                    if line[0] == '[':
                        break
                    hardstr += line + '  '
                    j += 1
            else:
                j = i + 1
            i = j

    return newstr, goodstr, hardstr



def func_save_data_file(words,str_=None,append=False):
    """save input words to data.txt without overwriting

    input:
        str OR 1D words: "word  word  ..." OR ["word", "word", ...]

        str_:       new, good, hard
        append:     whether save to a new file or backup

    output:
        format:
            [str_]
            word  word  word ...

    return:
        str, fname:  the name for backup
    OR,
        None:   when errors happen
    """
    wordslist = words.split() if isinstance(words,str) else words
    if len(wordslist) == 0:
        print('Warning: no inputs')
        return
    
    # outlines
    # format: ready to write to file
    if str_ is None: str_ = 'new'
    str_ = str_.lower()
    outlines = '[' + str_ + ']\n'

    line = ''
    for i in wordslist:
        line += i + '  '
        if len(line) >= 82:
            outlines += line.strip() + '\n'
            line = ''
    outlines += line.strip() + '\n\n\n'

    if not os.path.isdir('data'): os.makedirs('data')

    os.chdir('data')

    fname = ''
    if not append:
        if os.path.isfile('data.txt'):
            i = 1
            while True:
                fname = 'old-data-' + str(i) + '.txt'
                if not os.path.isfile(fname):
                    break
                i += 1
            # rename without any overwriting
            os.rename('data.txt',fname)

    # append writing
    with open('data.txt','a+') as f: f.write(outlines)

    os.chdir('../')

    return fname



class MLabel(tk.Label):
    def __init__(self,root,width,height,text=None,relief=None,bg=None,fg=None):
        self.image = tk.PhotoImage(width=width,height=height)
        relief = 'raised' if relief is None else relief
        tk.Label.__init__(self,root,text=text,relief=relief,
            image=self.image,compound='center',fg=fg,bg=bg
        )

class MButton(tk.Button):
    def __init__(self,root,width,height,text=None,relief=None,bg=None,fg=None,command=None):
        # it is important make img exist all the time
        self.image = tk.PhotoImage(width=width,height=height)
        relief = 'raised' if relief is None else relief
        tk.Button.__init__(self,root,text=text,relief=relief,
            image=self.image,compound='center',fg=fg,bg=bg,
            command=command
        )


class GUI(tk.Frame):
    def __init__(self,newlist=[],hardlist=[],goodlist=[]):
        self.root = tk.Tk()
        default_font = font.nametofont("TkDefaultFont")
        default_font.configure(size=16)

        self.newlist = newlist
        self.hardlist = hardlist
        self.goodlist = goodlist

        # for deletes
        self.deletelist = []

        # for finding previous
        self.pwordlist = []
        self.pindex = 0
        self.pflag = False

        self.index = 0
        self.reflist = []

        self.count = 0
        self.cur_button = None

        self.LWord = MLabel(self.root,600,250,text='Start')

        # function button
        self.BPrev = MButton(self.root,200,100,text='Prev',command=self.on_prev)
        BNext = MButton(self.root,200,100,text='Next',command=self.on_next)
        BDelete = MButton(self.root,200,50,text='Delete',command=self.on_delete)

        self.LInfo = MLabel(self.root,400,50,relief='flat',fg='indianred',text=self.count)
        self.LInfo['anchor'] = 'w'
        BReturn = MButton(self.root,200,50,text='Return',command=self.on_return)
        BSave = MButton(self.root,200,50,text='Save',command=self.on_save)

        self.BNew = MButton(self.root,200,100,text='New',command=self.on_new)
        self.BHard = MButton(self.root,200,100,text='Hard',command=self.on_hard)
        self.BGood = MButton(self.root,200,100,text='Good',command=self.on_good)
        LEntry = MLabel(self.root,200,100)
        SLabel = tk.Label(LEntry,text='Enter new word:')
        self.SEntry = tk.Entry(LEntry)
        SLabel.grid(row=0,column=0)
        self.SEntry.grid(row=1,column=0)
        self.SEntry.bind('<Return>',lambda t:self.on_entry())


        # align function button
        self.LWord.grid(row=0,rowspan=5,column=0,columnspan=3)
        BDelete.grid(row=0,column=3)
        self.BPrev.grid(row=1,rowspan=2,column=3)
        BNext.grid(row=3,rowspan=2,column=3)

        self.LInfo.grid(row=5,column=0,columnspan=2)
        BReturn.grid(row=5,column=2)
        BSave.grid(row=5,column=3)

        self.BNew.grid(row=6,column=0)
        self.BGood.grid(row=6,column=1)
        self.BHard.grid(row=6,column=2)
        LEntry.grid(row=6,column=3)

        self.root.protocol("WM_DELETE_WINDOW", self.on_exiting)
        self.root.mainloop()

    def on_exiting(self):
        if messagebox.askokcancel('Exiting', 'Don\'t forget to save!'):
            self.root.destroy()

    def on_save(self):
        """Save to data.txt file
        
        Be care: they may have repeats
        """
        bool_file = False
        if len(self.newlist) != 0:
            wlist = set(self.newlist)
            fname = func_save_data_file(wlist,'new')
            if not bool_file: bool_file = True if len(fname) != 0 else False

        if len(self.goodlist) != 0:
            wlist = set(self.goodlist)
            if bool_file:
                fname = func_save_data_file(wlist,'good',append=True)
            else:
                fname = func_save_data_file(wlist,'good')
            if not bool_file: bool_file = True if len(fname) != 0 else False

        if len(self.hardlist) != 0:
            wlist = set(self.hardlist)
            if bool_file:
                fname = func_save_data_file(wlist,'hard',append=True)
            else:
                fname = func_save_data_file(wlist,'hard')
            if not bool_file: bool_file = True if len(fname) != 0 else False

        if len(self.deletelist) != 0:
            wlist = set(self.deletelist)
            if bool_file:
                fname = func_save_data_file(wlist,'delete',append=True)
            else:
                fname = func_save_data_file(wlist,'delete')
            if not bool_file: bool_file = True if len(fname) != 0 else False

        info = ''
        if len(fname) != 0: info = 'Backup old file to ' + fname
        self.LInfo['text'] = 'Saved! ' + info


    def on_delete(self):
        """Delete the showing word, then move to the next"""
        if len(self.reflist) == 0:
            self.LInfo['text'] = 'No word is deleted!'
            return
        
        ndx = self.reflist[self.index]
        self.reflist.remove(ndx)

        if self.cur_button == 'new':
            wlist = self.newlist
        elif self.cur_button == 'good':
            wlist = self.goodlist
        elif self.cur_button == 'hard':
            wlist = self.hardlist
        else:
            return

        word = wlist[ndx]
        wlist.remove(word)
        if len(self.reflist) == 0:
            self.LWord['text'] = ''
            return
        
        self.reflist = [i-1 if i>ndx else i for i in self.reflist]
        self.LWord['text'] = wlist[self.reflist[self.index]]
        self.LInfo['text'] = 'Deleted word : ' + word
        self.deletelist.append(word)


    def on_entry(self):
        """Entry for new word"""
        word = self.SEntry.get()
        if len(word.split()) != 1:
            self.LInfo['text'] = 'Invalid new word!'
            return
        self.newlist.append(word)
        self.LInfo['text'] = 'Added new word: ' + word
        self.SEntry.delete(0,'end')

    def on_new(self):
        """Button New is pressed"""
        self.on_words_button('new')


    def on_good(self):
        """Button Good is pressed"""
        self.on_words_button('good')


    def on_hard(self):
        """Button Hard is pressed"""
        self.on_words_button('hard')


    def on_words_button(self,key):
        """Word Button is pressed, main function"""
        self.BPrev['state'] = 'normal'
        if key == 'new':
            btn = self.BNew
            wlist = self.newlist
        elif key == 'good':
            btn = self.BGood
            wlist = self.goodlist
        elif key == 'hard':
            btn = self.BHard
            wlist = self.hardlist
        else:
            return

        if self.cur_button is None:
            if key == 'new':
                self.BHard['text'] = 'Set to Hard'
                self.BGood['text'] = 'Set to Good'
            elif key == 'good':
                self.BHard['text'] = 'Set to Hard'
                self.BNew['text'] = 'Set to New'
            elif key == 'hard':
                self.BNew['text'] = 'Set to New'
                self.BGood['text'] = 'Set to Good'
            else:
                return

            btn['state'] = 'disabled'
            self.cur_button = key
            if len(wlist) == 0:
                self.reflist = []
                self.LInfo['text'] = 'No words!'
                return

            self.reflist = list(range(len(wlist)))
            random.shuffle(self.reflist)
            self.index = 0
            self.count += 1
            self.LWord['text'] = wlist[self.reflist[self.index]]
            self.LInfo['text'] = self.count
            self.pwordlist.append(self.LWord['text'])
            self.pflag = False
            self.wadd = 0
        else:
            if len(self.reflist) != 0:
                wlist.append(self.LWord['text'])
                self.wadd += 1


    def on_next(self):
        """Button Next is pressed"""
        self.BPrev['state'] = 'normal'
        if len(self.reflist) == 0:
            return
        
        self.index += 1
        self.count += 1
        self.LInfo['text'] = self.count
        if self.index >= len(self.reflist):
            self.index = 0
            random.shuffle(self.reflist)

        if self.cur_button == 'new':
            key = 'new'
            self.LWord['text'] = self.newlist[self.reflist[self.index]]
        elif self.cur_button == 'hard':
            key = 'hard'
            self.LWord['text'] = self.hardlist[self.reflist[self.index]]
        elif self.cur_button == 'good':
            key = 'good'
            self.LWord['text'] = self.goodlist[self.reflist[self.index]]
        else:
            return

        self.pwordlist.append(self.LWord['text'])
        self.pflag = False


    def on_prev(self):
        """Button Prev is pressed"""
        if not self.pflag:
            self.pindex = len(self.pwordlist) - 1
            self.pflag = True

        self.pindex -= 1
        if self.pindex < 0:
            self.BPrev['state'] = 'disabled'
            return
        self.BPrev['state'] = 'normal'
        self.LWord['text'] = self.pwordlist[self.pindex]
        self.LInfo['text'] = 'Reciting'
        
        


    def on_return(self):
        """Button Return is pressed"""
        if self.cur_button == 'new':
            self.BNew['state'] = 'normal'
        elif self.cur_button == 'good':
            self.BGood['state'] = 'normal'
        elif self.cur_button == 'hard':
            self.BHard['state'] = 'normal'
        
        self.BNew['text'] = 'New'
        self.BHard['text'] = 'Hard'
        self.BGood['text'] = 'Good'
        self.LInfo['text'] = self.count
        self.cur_button = None
        self.reflist = []




if __name__ == '__main__':

    # read newly added vocabulary.txt file
    newstr, hardstr, fnamelist = func_read_new_vocabulary()

    # save vocabulary.txt file without overwriting
    prodir = 'proceeded'
    if not os.path.isdir(prodir): os.makedirs(prodir)
    for file in fnamelist:
        cnt = 1
        while True:
            fnew = prodir + '/' + 'backup-' + str(cnt) + '-' + file
            if not os.path.isfile(fnew):
                break
            cnt += 1
        os.rename(file,fnew)

    newlist = newstr.split()
    hardlist = hardstr.split()
    if os.path.isdir('data'):
        os.chdir('data')
        newstr, goodstr, hardstr = func_read_data_file()
        newlist += newstr.split()
        hardlist += hardstr.split()
        goodlist = goodstr.split()
        os.chdir('../')
    
    GUI(newlist=newlist,hardlist=hardlist,goodlist=goodlist)

