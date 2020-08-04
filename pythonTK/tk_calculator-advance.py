import math
import tkinter as tk



def show(string):
    global g_exp, g_label, g_eval, g_tmp

    g_exp += string
    g_label['text'] = g_exp

    # replace human-readable exponential char ^ to computer-readable char **
    string = string.replace('^','**')
    g_tmp.append(string)
    
    if string == 'Ln(':
        g_eval += 'math.log('
    elif string in ('sin(','asin(','cos(','acos(','tan(','atan(',
                    'sqrt(','log10(','e**'):
        g_eval += 'math.' + string
    else:
        g_eval += string


def clear():
    global g_exp, g_label, g_rst, g_eval
    g_exp = ''
    g_label['text'] = g_exp
    g_eval = ''
    g_rst['text'] = ''


def deg():
    global g_deg, _pi
    if g_deg.get() == 'Deg':
        g_deg.set('Rad')
        _pi = 1.0
    else:
        g_deg.set('Deg')
        _pi = 180.0


def fdel():
    global g_exp, g_label, g_eval, g_tmp

    if len(g_tmp) != 0:

        # delete the evaluation expression
        string = g_tmp.pop()
        if string == 'Ln(':
            g_eval = g_eval[:g_eval.rfind('math.log(')]
        elif string in ('sin(','asin(','cos(','acos','tan(',
                        'atan(','sqrt(','log10(','e**'):
            g_eval = g_eval[:g_eval.rfind('math.' + string)]
        else:
            g_eval = g_eval[:g_eval.rfind(string)]

        # delete the showing expression
        # replace computer-readable char ** to human-readable exponential char ^
        tmp = string.replace('**','^')
        g_exp = g_exp[:g_exp.rfind(tmp)]
        g_label['text'] = g_exp


           
def calc():
    global g_eval, g_rst, _pi
  
    def proeval(sub,string):
        # separately evaluate trigonometric function
        ref = string.find(sub)
        if ref != -1:
            start = string[ref:].find('(')
            if ref + start >= len(string): return 'ERR'
            
            end = string[ref+start:].find(')')
            if end != -1:
                try:                
                    rst = eval(string[ref+start+1:ref+start+end])                
                    rst = rst / _pi * math.pi

                    rst = string[ref:ref+start+1] + str(rst) + ')'
                    rst = eval(rst)

                    if ref + start + end + 1 >= len(string):
                        t = ''
                    else:
                        t = string[ref+start+end+1:]
                        
                    string = string[:ref] + '(' + str(rst) + ')' + t
                except:
                    string = 'ERR'
            else:
                string = 'ERR'
        return string
    

    bo = False
    result = g_eval
    for tri in ('sin(','asin(','cos(','acos(','tan(','atan('):
        while True:
            if result.find('math.'+tri) == -1:
                break
            else:
                result = proeval('math.'+tri,result)
                if result == 'ERR':
                    bo = True
                    break
        if bo: break
    if bo:
        g_rst['text'] = 'ERROR'
    else:
        try:
            result = eval(result)
            g_rst['text'] = result
        except:
            g_rst['text'] = 'ERROR'


# global windows parameter
g_width = 320
g_height = 200
g_title = 'PythonicTk Scientific Calculator'

# global, expression
g_exp = ''

# global, result
g_eval = ''
_pi = 180

# global, used for delete/backspace
g_tmp = []



root = tk.Tk()
root.geometry('{:}x{:}'.format(g_width,g_height))
#root.resizable(0,0)
root.title(g_title)


# A trick to set tk.Button in pixel
t_width = int(g_width / 6) - 10
t_height = int((g_height-10) / 8) - 10

g_pixel = tk.PhotoImage(width=t_width,height=t_height)


# Window to show expression
frame_exp = tk.Frame(root,width=g_width,height=t_height+5)
frame_exp.grid(row=0,column=0,columnspan=6)
g_label = tk.Label(frame_exp,text='',width=g_width,anchor=tk.E)
g_label.place(width=g_width,height=t_height+5)


# Window to show result
frame_rst = tk.Frame(root,width=g_width,height=t_height+5,bd=1,highlightbackground='black',
                     highlightcolor='black',highlightthickness=1)
frame_rst.grid(row=1,column=0,columnspan=6)
g_rst = tk.Label(frame_rst,text='',width=g_width,anchor=tk.E)
g_rst.place(width=g_width,height=t_height+5)

# global, control degree or radius
g_deg = tk.StringVar()
g_deg.set('Deg')



# Layout

# sqrt, %, 10^n, Log10, clear, del
tk.Button(root,image=g_pixel,compound=tk.CENTER,bd=3,text='sqrt(',command=lambda:show('sqrt(')).grid(row=2,column=0)
tk.Button(root,image=g_pixel,compound=tk.CENTER,bd=3,text='%',command=lambda:show('%')).grid(row=2,column=1)
tk.Button(root,image=g_pixel,compound=tk.CENTER,bd=3,text='10^',command=lambda:show('10^')).grid(row=2,column=2)
tk.Button(root,image=g_pixel,compound=tk.CENTER,bd=3,text='Log10(',command=lambda:show('log10(')).grid(row=2,column=3)
tk.Button(root,image=g_pixel,compound=tk.CENTER,bd=3,text='clear',command=lambda:clear()).grid(row=2,column=4)
tk.Button(root,image=g_pixel,compound=tk.CENTER,bd=3,text='del',command=lambda:fdel()).grid(row=2,column=5)


# ^2, ^, e^n, Ln, sin, asin
tk.Button(root,image=g_pixel,compound=tk.CENTER,bd=3,text='^2',command=lambda:show('^2')).grid(row=3,column=0)
tk.Button(root,image=g_pixel,compound=tk.CENTER,bd=3,text='^',command=lambda:show('^')).grid(row=3,column=1)
tk.Button(root,image=g_pixel,compound=tk.CENTER,bd=3,text='e^',command=lambda:show('e^')).grid(row=3,column=2)
tk.Button(root,image=g_pixel,compound=tk.CENTER,bd=3,text='Ln(',command=lambda:show('Ln(')).grid(row=3,column=3)
tk.Button(root,image=g_pixel,compound=tk.CENTER,bd=3,text='sin(',command=lambda:show('sin(')).grid(row=3,column=4)
tk.Button(root,image=g_pixel,compound=tk.CENTER,bd=3,text='asin(',command=lambda:show('asin(')).grid(row=3,column=5)


# 1, 2, 3, /, cos, acos
tk.Button(root,image=g_pixel,compound=tk.CENTER,bd=3,text='1',command=lambda:show('1')).grid(row=4,column=0)
tk.Button(root,image=g_pixel,compound=tk.CENTER,bd=3,text='2',command=lambda:show('2')).grid(row=4,column=1)
tk.Button(root,image=g_pixel,compound=tk.CENTER,bd=3,text='3',command=lambda:show('3')).grid(row=4,column=2)
tk.Button(root,image=g_pixel,compound=tk.CENTER,bd=3,text='/',command=lambda:show('/')).grid(row=4,column=3)
tk.Button(root,image=g_pixel,compound=tk.CENTER,bd=3,text='cos(',command=lambda:show('cos(')).grid(row=4,column=4)
tk.Button(root,image=g_pixel,compound=tk.CENTER,bd=3,text='acos(',command=lambda:show('acos(')).grid(row=4,column=5)


# 4, 5, 6, *, tan, atan
tk.Button(root,image=g_pixel,compound=tk.CENTER,bd=3,text='4',command=lambda:show('4')).grid(row=5,column=0)
tk.Button(root,image=g_pixel,compound=tk.CENTER,bd=3,text='5',command=lambda:show('5')).grid(row=5,column=1)
tk.Button(root,image=g_pixel,compound=tk.CENTER,bd=3,text='6',command=lambda:show('6')).grid(row=5,column=2)
tk.Button(root,image=g_pixel,compound=tk.CENTER,bd=3,text='*',command=lambda:show('*')).grid(row=5,column=3)
tk.Button(root,image=g_pixel,compound=tk.CENTER,bd=3,text='tan(',command=lambda:show('tan(')).grid(row=5,column=4)
tk.Button(root,image=g_pixel,compound=tk.CENTER,bd=3,text='atan(',command=lambda:show('atan(')).grid(row=5,column=5)



# 7, 8, 9, +, (, ),
tk.Button(root,image=g_pixel,compound=tk.CENTER,bd=3,text='7',command=lambda:show('7')).grid(row=6,column=0)
tk.Button(root,image=g_pixel,compound=tk.CENTER,bd=3,text='8',command=lambda:show('8')).grid(row=6,column=1)
tk.Button(root,image=g_pixel,compound=tk.CENTER,bd=3,text='9',command=lambda:show('9')).grid(row=6,column=2)
tk.Button(root,image=g_pixel,compound=tk.CENTER,bd=3,text='+',command=lambda:show('+')).grid(row=6,column=3)
tk.Button(root,image=g_pixel,compound=tk.CENTER,bd=3,text='(',command=lambda:show('(')).grid(row=6,column=4)
tk.Button(root,image=g_pixel,compound=tk.CENTER,bd=3,text=')',command=lambda:show(')')).grid(row=6,column=5)


# 0, ., E, -, rad/deg, =
tk.Button(root,image=g_pixel,compound=tk.CENTER,bd=3,text='0',command=lambda:show('0')).grid(row=7,column=0)
tk.Button(root,image=g_pixel,compound=tk.CENTER,bd=3,text='.',command=lambda:show('.')).grid(row=7,column=1)
tk.Button(root,image=g_pixel,compound=tk.CENTER,bd=3,text='E',command=lambda:show('E')).grid(row=7,column=2)
tk.Button(root,image=g_pixel,compound=tk.CENTER,bd=3,text='-',command=lambda:show('-')).grid(row=7,column=3)
tk.Button(root,image=g_pixel,compound=tk.CENTER,bd=3,textvariable=g_deg,command=lambda:deg()).grid(row=7,column=4)
tk.Button(root,image=g_pixel,compound=tk.CENTER,bd=3,text='=',command=lambda:calc()).grid(row=7,column=5)

root.mainloop()
