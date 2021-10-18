import sys
import PySimpleGUI as sg
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk


def func(x, y):
    return 3 * y ** (2 / 3)
#     return (y ** 2 + x * y - x ** 2) / (x ** 2)
#     return (y / x + x * np.cos(x))


def exact_function(x):
    return (x - 1) ** 3
#     return x * (1 + (x ** 2) / 3) / (1 - (x ** 2) / 3)
#     return (x * np.sin(x) + x / np.pi)


class ExactSolution:
    @staticmethod
    def Exact(h, y0, x0, X, N):
        x = x0
        y_vals = [exact_function(x)]
        
        for i in range(N):
            x = x + h
            y_vals.append(exact_function(x))
        
        return y_vals

class GlobalErrors:
    @staticmethod
    def Global_Errors(exact_values, given_values):
        e_vals = []
        for exact, given in zip(exact_values, given_values):
            e_vals.append(abs(exact - given))
        return e_vals


class AbstractMethod:
    def inner_function(self, x_cur, y_cur): raise NotImlementedError

    def Method(self): raise NotImplementedError

    def Work_Method(self): raise NotImplementedError

    def All_Method(self): raise NotImplementedError

    def LocarErrors(self, y_exact, x_values): raise NotImplementedError


class EulerMethod(AbstractMethod):
    def __init__(self, h, x0, y0, X, N):
        self.h = h
        self.x0 = x0
        self.y0 = y0
        self.X = X
        self.N = N
    
    
    def inner_function(self, x_cur, y_cur):
        return self.h * func(x_cur, y_cur)
    

    def Method(self, x_values):
        y_vals = [self.y0]
        y_cur = self.y0
        
        for i in range(1, len(x_values)):
            y_cur = y_cur + self.inner_function(x_values[i], y_cur)
            y_vals.append(y_cur)
            
        return y_vals
        
        
    def Local_Errors(self, y_exact, x_values):
        e_vals = [0.0]
        
        for i in range(1, len(x_values)):
            e_vals.append(abs(y_exact[i] - y_exact[i - 1] -\
                                self.inner_function(x_cur=x_values[i - 1], y_cur=y_exact[i - 1])))
        return e_vals


    def Work_Method(self):
        x_values = [self.x0]
        
        for i in range(self.N):
            x_values.append(x_values[-1] + self.h)

        euler_values = self.Method(x_values)

        exact_values = ExactSolution.Exact(self.h, self.y0, self.x0, self.X, self.N)

        global_error_values = GlobalErrors.Global_Errors(exact_values, euler_values)

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
        increasment = self.h * (k1 + k2) / 2
        return increasment
    
    
    def Method(self, x_values):
        y_vals = [self.y0]
        y_cur = self.y0
        
        for i in range(1, len(x_values)):
            y_cur = y_cur + self.inner_function(x_values[i], y_cur)
            y_vals.append(y_cur)
            
        return y_vals
    
    
    def Local_Errors(self, y_exact, x_values):
        e_vals = [0.0]
        
        for i in range(1, len(x_values)):
            e_vals.append(abs(y_exact[i] - y_exact[i - 1] -\
                            self.inner_function(x_cur = x_values[i - 1], y_cur = y_exact[i - 1])))
        return e_vals
    
    
    def Work_Method(self):
        x_values = [self.x0]
        
        for i in range(self.N):
            x_values.append(x_values[-1] + self.h)

        improved_euler_values = self.Method(x_values)

        exact_values = ExactSolution.Exact(self.h, self.y0, self.x0, self.X, self.N)

        global_error_values = GlobalErrors.Global_Errors(exact_values, improved_euler_values)

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
        k1 = func(x = x_cur,                y = y_cur)
        k2 = func(x = (x_cur + self.h / 2), y = (y_cur + self.h * k1 / 2))
        k3 = func(x = (x_cur + self.h / 2), y = (y_cur + self.h * k2 / 2))
        k4 = func(x = (x_cur + self.h),     y = (y_cur + self.h * k3))

        increasment = self.h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
        return increasment
        
        
    def Method(self, x_values):
        y_vals = [self.y0]
        y_cur = self.y0
        
        for i in range(1, len(x_values)):
            y_cur = y_cur + self.inner_function(x_values[i], y_cur)
            y_vals.append(y_cur)
            
        return y_vals
    
    
    def Local_Errors(self, y_exact, x_values):
        e_vals = [0.0]
        
        for i in range(1, len(x_values)):
            e_vals.append(abs(y_exact[i] - y_exact[i - 1] -\
                        self.inner_function(x_cur = x_values[i - 1], y_cur = y_exact[i - 1])))
        return e_vals

        
    def Work_Method(self):
        x_values = [self.x0]
        
        for i in range(self.N):
            x_values.append(x_values[-1] + self.h)

        runge_kutta_values = self.Method(x_values)

        exact_values = ExactSolution.Exact(self.h, self.y0, self.x0, self.X, self.N)

        global_error_values = GlobalErrors.Global_Errors(exact_values, runge_kutta_values)

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
        [sg.Text('Choose task:'), sg.InputCombo(['Approximation', 'LTE', 'GTE', 'GTE for N'],\
                                                  default_value = 'Approximation', size = (20, 4))],
        [sg.Text('Mark methods to use:'), sg.Checkbox('Euler'), sg.Checkbox('Improved Euler'),\
            sg.Checkbox('Runge-Kutta')],
        [sg.Text('x0', size = (2, 1)), sg.InputText(size = (10, 1))],
        [sg.Text('y0', size = (2, 1)), sg.InputText(size = (10, 1))],
        [sg.Text('X',  size = (2, 1)), sg.InputText(size = (10, 1))],
        [sg.Text('N',  size = (2, 1)), sg.InputText(size = (10, 1))],
        [sg.Button('Plot')],
        [sg.Canvas(key = 'controls_cv')],
        [sg.Text('Figure:')],
        [sg.Column(
            layout = [
                [sg.Canvas(key = 'fig_cv',
                           # size of the gui
                           size = (400 * 2, 300)
                           )]
            ],
            background_color = '#DAE0E6', 
            pad = (0, 0)
        )],
        [sg.Button('Exit')]
    ]

    window = sg.Window('DE approximation with controls', layout) # '#e7f4f4'

    while True:
        event, values = window.read()
        
        if event in (sg.WIN_CLOSED, 'Exit'):
            break
            
        is_error = 0
        
        task = values[0]
        is_Euler = values[1]
        is_Improved_Euler = values[2]
        is_RungeKutta = values[3]
        x0 = str(values[4])
        
        try:
            x0 = float(x0.replace(',', '.'))
        except:
            if not is_error: 
                sg.popup('All values should be float or integer, N > 0\nChange x0')
                is_error = 1

        y0 = str(values[5])
        
        try:
            y0 = float(y0.replace(',', '.'))
        except:
            if not is_error: 
                sg.popup('All values should be float or integer, N > 0\nChange y0')
                is_error = 1
            
        X  = str(values[6])
        
        try:
            X = float(X.replace(',', '.'))
        except:
            if not is_error: 
                sg.popup('All values should be float or integer, N > 0\nChange X')
                is_error = 1
        
        N = 1
        
        try:
            N = int(values[7])
        except:
            if not is_error: 
                sg.popup('N should be greater than zero and integer')
                is_error = 1
        
        if N == 0:
            if not is_error: 
                sg.popup('N should be greater than zero and integer')
                is_error = 1
            
        if x0 == X:
            if not is_error: 
                sg.popup('X must NOT be equal to x0')
                is_error = 1
        
        if is_error:
            continue
        
        h  = X - x0
        
        try:
            h = float((X - x0) / N)
        except:
            h = X - x0
        
        euler          = EulerMethod(h, x0, y0, X, N)
        improved_euler = ImprovedEulerMethod(h, x0, y0, X, N)
        runge_kutta    = RungeKuttaMethod(h, x0, y0, X, N)
            
        euler_values          = euler.All_Method()
        improved_euler_values = improved_euler.All_Method()
        runge_kutta_values    = runge_kutta.All_Method()
        
#       return (x_values, runge_kutta_values, exact_values, local_error_values, global_error_values)

        if event == 'Plot':
            check_local_error = 0
            
            if task == 'GTE for N':
                if not is_Euler and not is_Improved_Euler and not is_RungeKutta:
                    sg.popup('Choose method for which you want to plot GTE')
                    check_local_error = 1
            
            if check_local_error:
                continue
            
            plt.clf()
            plt.figure(1)
            fig = plt.gcf()
            DPI = fig.get_dpi()

            title = task
            x_label = 'X'
            y_label = 'Y'
        
            if task == 'Approximation':
                x_exact = np.linspace(x0, X, 1000)
                y_exact = np.power(np.subtract(x_exact, 1), 3)

                plt.plot(x_exact, y_exact, 'g', label='exact solution')

                if is_Euler:
                    plt.plot(euler_values[0], euler_values[1], 'b',label='Euler')

                if is_Improved_Euler:
                    plt.plot(improved_euler_values[0], improved_euler_values[1], 'y',label='Improved Euler')

                if is_RungeKutta:
                    plt.plot(runge_kutta_values[0], runge_kutta_values[1], 'r',label='Runge Kutta')
            
            elif task == 'LTE':
                if is_Euler:
                    plt.plot(euler_values[0], euler_values[3], 'b',label='Euler LTE')

                if is_Improved_Euler:
                    plt.plot(improved_euler_values[0], improved_euler_values[3], 'y',label='Improved Euler LTE')

                if is_RungeKutta:
                    plt.plot(runge_kutta_values[0], runge_kutta_values[3], 'r',label='Runge Kutta LTE')
    
            
            elif task == 'GTE':
                if is_Euler:
                    plt.plot(euler_values[0], euler_values[4], 'b',label='Euler GTE')

                if is_Improved_Euler:
                    plt.plot(improved_euler_values[0], improved_euler_values[4], 'y',label='Improved Euler GTE')

                if is_RungeKutta:
                    plt.plot(runge_kutta_values[0], runge_kutta_values[4], 'r',label='Runge Kutta GTE')
                    
            elif task == 'GTE for N':
                title = 'GTE for N = [1, 100]'
                
                euler_GTE = []
                improved_euler_GTE = []
                runge_kutta_GTE = []
                
                N_values = list(range(1, 101))
                for i in N_values:
                    h_2 = 1
                    try:
                        h_2 = (X - x0) / i
                    except ZeroDivisionError:
                        pass
                    
                    euler_2          = EulerMethod(h_2, x0, y0, X, i)
                    improved_euler_2 = ImprovedEulerMethod(h_2, x0, y0, X, i)
                    runge_kutta_2    = RungeKuttaMethod(h_2, x0, y0, X, i)
            
                    euler_values_2          = euler_2.All_Method()
                    improved_euler_values_2 = improved_euler_2.All_Method()
                    runge_kutta_values_2    = runge_kutta_2.All_Method()
                    
                    euler_GTE.append(max(euler_values_2[4]))
                    improved_euler_GTE.append(max(improved_euler_values_2[4]))
                    runge_kutta_GTE.append(max(runge_kutta_values_2[4]))
                
                if is_Euler:
                    plt.plot(N_values, euler_GTE, 'b',label='Euler GTE')

                if is_Improved_Euler:
                    plt.plot(N_values, improved_euler_GTE, 'y',label='Improved Euler GTE')

                if is_RungeKutta:
                    plt.plot(N_values, runge_kutta_GTE, 'r',label='Runge Kutta GTE')
                    
                x_label = 'N'
                y_label = 'Max GTE'
                    
            if check_local_error:
                continue
                
            plt.title(f'{title}')
            plt.xlabel(x_label)
            plt.ylabel(y_label)
            plt.grid()
            plt.legend(fontsize = 15, loc = 0)

            draw_figure_w_toolbar(window['fig_cv'].TKCanvas, fig, window['controls_cv'].TKCanvas)
        
        else:
            sg.Popup('Ok clicked', keep_on_top = True)

    window.close()


if __name__ == '__main__':
    orig_stdout = sys.stdout
    f = open('output.txt', 'w')
    sys.stdout = f

    gui_trial()

    sys.stdout = orig_stdout
    f.close()
