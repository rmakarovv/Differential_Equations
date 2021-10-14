import sys
import PySimpleGUI as sg
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

Rounding_Number = 10

def func(x, y):
#     return (y ** 2 + x * y - x ** 2) / (x ** 2)
    return 3 * y ** (2 / 3)


def exact_function(x):
#     return x * (1 + (x ** 2) / 3) / (1 - (x ** 2) / 3)
    return (x - 1) ** 3


class ExactSolution:
    @staticmethod
    def Exact(h, y0, x0, X, N):
        x = x0
        y_vals = [exact_function(x)]
        
        for i in range(N):
            x = round(x + h, Rounding_Number)
            y_vals.append(exact_function(x))
        
        return y_vals


def Global_Errors(exact_values, given_values):
    e_vals = []
    for exact, given in zip(exact_values, given_values):
        e_vals.append(round(exact - given, Rounding_Number))
    return e_vals


class AbstractMethod:
    def Method(self): raise NotImplementedError

    def Work_Method(self): raise NotImplementedError

    def All_Method(self): raise NotImplementedError
        
    def inner_function(self, x_cur, y_cur): raise NotImlementedError
    
    def LocarErrors(self, y_exact, x_values): raise NotImplementedError


class EulerMethod(AbstractMethod):
    def __init__(self, h, x0, y0, X, N):
        self.h = h
        self.x0 = x0
        self.y0 = y0
        self.X = X
        self.N = N
    
    
    def inner_function(self, x_cur, y_cur):
        return round(y_cur + self.h * func(x_cur, y_cur), Rounding_Number)
    

    def Method(self, x_values):
        y_vals = [self.y0]
        y_cur = self.y0
        
        for i in range(1, len(x_values)):
            y_cur = self.inner_function(x_values[i], y_cur)
            y_vals.append(y_cur)
            
        return y_vals
        
        
    def Local_Errors(self, y_exact, x_values):
        e_vals = [0.0]
        
        for i in range(1, len(x_values)):
            e_vals.append(round(y_exact[i] - y_exact[i - 1] - self.h * \
                                self.inner_function(x_cur=x_values[i - 1], y_cur=y_exact[i - 1]), Rounding_Number))
        return e_vals


    def Work_Method(self):
        x_values = [self.x0]
        
        for i in range(self.N):
            x_values.append(round(x_values[-1] + self.h, Rounding_Number))

        euler_values = self.Method(x_values)

        exact_values = ExactSolution.Exact(self.h, self.y0, self.x0, self.X, self.N)

        global_error_values = Global_Errors(exact_values, euler_values)

        local_error_values = self.Local_Errors(exact_values, x_values)

        return (x_values, euler_values, exact_values, local_error_values, global_error_values)


    def All_Method(self):
        thing = self.Work_Method()
        return thing



class ImprovedEulerMethod(AbstractMethod):
    def __init__(self, h, x0, y0, X, N):
        self.h = h
        self.x0 = x0
        self.y0 = y0
        self.X = X
        self.N = N

        
    def inner_function(self, x_cur, y_cur):
        k1 = func(x_cur, y_cur)
        k2 = func(x_cur + self.h, y_cur + self.h * k1)
        y = round(y_cur + self.h * (k1 + k2) / 2, Rounding_Number)
        return y
    
    
    def Method(self, x_values):
        y_cur = self.y0
        y_vals = [self.y0]
        
        for i in range(1, len(x_values)):
            y_cur = self.inner_function(x_cur = x_values[i], y_cur = y_cur)
            y_vals.append(y_cur)

        return y_vals
    
    
    def Local_Errors(self, y_exact, x_values):
        e_vals = [0.0]
        
        for i in range(1, len(x_values)):
            e_vals.append(round(y_exact[i] - y_exact[i - 1] - self.h * \
                                self.inner_function(x_cur = x_values[i - 1], y_cur = y_exact[i - 1]), Rounding_Number))
        return e_vals
    
    
    def Work_Method(self):
        x_values = [self.x0]
        
        for i in range(self.N):
            x_values.append(round(x_values[-1] + self.h, Rounding_Number))

        improved_euler_values = self.Method(x_values)

        exact_values = ExactSolution.Exact(self.h, self.y0, self.x0, self.X, self.N)

        global_error_values = Global_Errors(exact_values, improved_euler_values)

        local_error_values = self.Local_Errors(exact_values, x_values)

        return (x_values, improved_euler_values, exact_values, local_error_values, global_error_values)

    
    def All_Method(self):
        values = self.Work_Method()
        return values


class RungeKuttaMethod(AbstractMethod):
    def __init__(self, h, x0, y0, X, N):
        self.h = h
        self.x0 = x0
        self.y0 = y0
        self.X = X
        self.N = N
        
    def inner_function(self, x_cur, y_cur):
        k1 = func(x=x_cur, y=y_cur)
        k2 = func(x=(x_cur + self.h / 2), y=(y_cur + self.h * k1 / 2))
        k3 = func(x=(x_cur + self.h / 2), y=(y_cur + self.h * k2 / 2))
        k4 = func(x=(x_cur + self.h), y=(y_cur + self.h * k3))

        y_cur = round(y_cur + self.h * (k1 + 2 * k2 + 2 * k3 + k4) / 6, Rounding_Number)
        return y_cur
        
        
    def Method(self, x_values):
        y_cur = self.y0
        y_vals = [self.y0]
        
        for i in range(1, len(x_values)):
            y_cur = self.inner_function(x_values[i], y_cur)
            y_vals.append(y_cur)

        return y_vals
    
    
    def Local_Errors(self, y_exact, x_values):
        e_vals = [0.0]
        
        for i in range(1, len(x_values)):
            e_vals.append(round(y_exact[i] - y_exact[i - 1] - self.h * \
                                self.inner_function(x_cur = x_values[i - 1], y_cur = y_exact[i - 1]), Rounding_Number))
        return e_vals

        
    def Work_Method(self):
        x_values = [self.x0]
        
        for i in range(self.N):
            x_values.append(round(x_values[-1] + self.h, Rounding_Number))

        runge_kutta_values = self.Method(x_values)

        exact_values = ExactSolution.Exact(self.h, self.y0, self.x0, self.X, self.N)

        global_error_values = Global_Errors(exact_values, runge_kutta_values)

        local_error_values = self.Local_Errors(exact_values, x_values)

        return (x_values, runge_kutta_values, exact_values, local_error_values, global_error_values)


    def All_Method(self):
        values = self.Work_Method()
        return values


def draw_figure_w_toolbar(canvas, fig, canvas_toolbar):
    if canvas.children:
        for child in canvas.winfo_children():
            child.destroy()
    if canvas_toolbar.children:
        for child in canvas_toolbar.winfo_children():
            child.destroy()
    figure_canvas_agg = FigureCanvasTkAgg(fig, master=canvas)
    figure_canvas_agg.draw()
    toolbar = Toolbar(figure_canvas_agg, canvas_toolbar)
    toolbar.update()
    figure_canvas_agg.get_tk_widget().pack(side='right', fill='both', expand=1)


class Toolbar(NavigationToolbar2Tk):
    def __init__(self, *args, **kwargs):
        super(Toolbar, self).__init__(*args, **kwargs)
        

        
def gui_trial():
    layout = [
        [sg.Text('Choose Method:'), sg.InputCombo(['Euler', 'Improved Euler', 'Runge-Kutta'],\
                                                  default_value ='Euler', size=(20, 3))],
        [sg.Text('x0'), sg.InputText(size=(10, 1))],
        [sg.Text('y0'), sg.InputText(size=(10, 1))],
        [sg.Text(' X'), sg.InputText(size=(10, 1))],
        [sg.Text(' N'), sg.InputText(size=(10, 1))],
        [sg.Button('Plot')],
        [sg.Canvas(key='controls_cv')],
        [sg.Text('Figure:')],
        [sg.Column(
            layout=[
                [sg.Canvas(key='fig_cv',
                           # it's important that you set this size
                           size=(400 * 2, 300)
                           )]
            ],
            background_color='#DAE0E6',
            pad=(0, 0)
        )],
        [sg.Button('Exit')]
    ]

    window = sg.Window('DE approximation with controls', layout)

    while True:
        event, values = window.read()
        
        if event in (sg.WIN_CLOSED, 'Exit'):
            break
        
        method = values[0]
        x0 = float(values[1])
        y0 = float(values[2])
        X = float(values[3])
        N = int(values[4])
        h = round((X - x0) / N, Rounding_Number)
        
        
        Method = EulerMethod(h, x0, y0, X, N)
        
        if method == 'Improved Euler':
            Method = ImprovedEulerMethod(h, x0, y0, X, N)
        elif method == 'Runge-Kutta':
            Method = RungeKuttaMethod(h, x0, y0, X, N)
            
        values = Method.All_Method()
        x_vals = values[0]
        euler_vals = values[1]
        exact_vals = values[2]
        local_error_vals = values[3]
        global_error_vals = values[4]

        if event == 'Plot':
            plt.clf()
            # ------------------------------- PASTE YOUR MATPLOTLIB CODE HERE
            plt.figure(3)
            fig = plt.gcf()
            DPI = fig.get_dpi()


#             # ------------------------------- you have to play with this size to reduce the 
#             movement error when the mouse hovers over the figure, it's close to canvas size
            fig.set_size_inches(404 * 2 / float(DPI), 404 / float(DPI))


            x = x_vals
            y = euler_vals
            y_exact = exact_vals

            plt.plot(x, y, 'b',label='approximation')
            plt.plot(x, y_exact, 'g', label='exact solution')
            
            plt.title(f'{method} method approximation')
            plt.xlabel('X')
            plt.ylabel('Y')
            plt.grid()
            plt.legend(fontsize=15, loc=0)
#             fig = plt.figure()

            # ------------------------------- Instead of plt.show()
            draw_figure_w_toolbar(window['fig_cv'].TKCanvas, fig, window['controls_cv'].TKCanvas)
        else:
            sg.Popup('Ok clicked', keep_on_top=True)

    window.close()


if __name__ == '__main__':
    orig_stdout = sys.stdout
    f = open('output.txt', 'w')
    sys.stdout = f

    gui_trial()

    sys.stdout = orig_stdout
    f.close()
