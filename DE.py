#
#  Numerical methods for solving DE and GUI with graph plotting
#  implementation.
#
#  @author Roman Makarov, BS20_06, o.makarov@innopolis.university
#

import sys
import PySimpleGUI as sg
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk


class Functions:
    """Interface for functions"""
    @staticmethod
    def function(x, y): raise NotImplementedError

    @staticmethod
    def exact_function(x): raise NotImplementedError


class MyFunctions(Functions):
    @staticmethod
    def function(x, y):
        return 3 * y ** (2 / 3)

    @staticmethod
    def exact_function(x):
        return (x - 1) ** 3


class ExactSolution:
    @staticmethod
    def exact(h, x0, N):
        x = x0
        y_vals = [MyFunctions.exact_function(x)]

        for i in range(N):
            x = x + h
            y_vals.append(MyFunctions.exact_function(x))

        return y_vals


class LocalErrors:
    @staticmethod
    def local_Errors(exact_values, given_values):
        e_vals = []
        for exact, given in zip(exact_values, given_values):
            e_vals.append(abs(exact - given))
        return e_vals


class NumericalMethod:
    """Interface for numerical methods"""
    def inner_function(self, x_cur, y_cur): raise NotImplementedError

    def method(self, x_values): raise NotImplementedError

    def work_method(self): raise NotImplementedError


class Euler(NumericalMethod):
    def __init__(self, h, x0, y0, X, N):
        self.h = h
        self.x0 = x0
        self.y0 = y0
        self.X = X
        self.N = N

    def inner_function(self, x_cur, y_cur):
        return self.h * MyFunctions.function(x=x_cur, y=y_cur)

    def method(self, x_values):
        y_vals = [self.y0]
        y_cur = self.y0

        for i in range(1, len(x_values)):
            y_cur = y_cur + self.inner_function(x_cur=x_values[i], y_cur=y_cur)
            y_vals.append(y_cur)

        return y_vals

    def work_method(self):
        x_values = [self.x0]

        for i in range(self.N):
            x_values.append(x_values[-1] + self.h)

        euler_values = self.method(x_values)

        exact_values = ExactSolution.exact(h=self.h, x0=self.x0, N=self.N)

        local_error_values = LocalErrors.local_Errors(exact_values, euler_values)

        global_error_value = max(local_error_values)

        return (x_values, euler_values, exact_values, local_error_values, global_error_value)


class ImprovedEuler(NumericalMethod):
    def __init__(self, h, x0, y0, X, N):
        self.h = h
        self.x0 = x0
        self.y0 = y0
        self.X = X
        self.N = N

    def inner_function(self, x_cur, y_cur):
        k1 = MyFunctions.function(x=x_cur, y=y_cur)
        k2 = MyFunctions.function(x=(x_cur + self.h), y=(y_cur + self.h * k1))
        increment = self.h * (k1 + k2) / 2
        return increment

    def method(self, x_values):
        y_vals = [self.y0]
        y_cur = self.y0

        for i in range(1, len(x_values)):
            y_cur = y_cur + self.inner_function(x_cur=x_values[i], y_cur=y_cur)
            y_vals.append(y_cur)

        return y_vals

    def work_method(self):
        x_values = [self.x0]

        for i in range(self.N):
            x_values.append(x_values[-1] + self.h)

        improved_euler_values = self.method(x_values)

        exact_values = ExactSolution.exact(h=self.h, x0=self.x0, N=self.N)

        local_error_values = LocalErrors.local_Errors(exact_values, improved_euler_values)

        global_error_value = max(local_error_values)

        return (x_values, improved_euler_values, exact_values, local_error_values, global_error_value)


class RungeKutta(NumericalMethod):
    def __init__(self, h, x0, y0, X, N):
        self.h = h
        self.x0 = x0
        self.y0 = y0
        self.X = X
        self.N = N

    def inner_function(self, x_cur, y_cur):
        k1 = MyFunctions.function(x=x_cur,                  y=y_cur)
        k2 = MyFunctions.function(x=(x_cur + self.h / 2),   y=(y_cur + self.h * k1 / 2))
        k3 = MyFunctions.function(x=(x_cur + self.h / 2),   y=(y_cur + self.h * k2 / 2))
        k4 = MyFunctions.function(x=(x_cur + self.h),       y=(y_cur + self.h * k3))

        increment = self.h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
        return increment

    def method(self, x_values):
        y_vals = [self.y0]
        y_cur = self.y0

        for i in range(1, len(x_values)):
            y_cur = y_cur + self.inner_function(x_values[i], y_cur)
            y_vals.append(y_cur)

        return y_vals

    def work_method(self):
        x_values = [self.x0]

        for i in range(self.N):
            x_values.append(x_values[-1] + self.h)

        runge_kutta_values = self.method(x_values)

        exact_values = ExactSolution.exact(h=self.h, x0=self.x0, N=self.N)

        local_error_values = LocalErrors.local_Errors(exact_values, runge_kutta_values)

        global_error_value = max(local_error_values)

        return (x_values, runge_kutta_values, exact_values, local_error_values, global_error_value)


# Toolbar for figure
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
    figure_canvas_agg.get_tk_widget().pack(side='top', fill='both', expand=1)


class Toolbar(NavigationToolbar2Tk):
    def __init__(self, *args, **kwargs):
        super(Toolbar, self).__init__(*args, **kwargs)


def gui():
    layout = [
        [sg.Text('Choose task:'), sg.InputCombo(['Approximation', 'LTE', 'GTE for N'], \
                                                default_value='Approximation', size=(20, 4))],
        [sg.Text('Mark methods to use:'), sg.Checkbox('Euler'), sg.Checkbox('Improved Euler'), \
         sg.Checkbox('Runge-Kutta')],
        [sg.Text('x0',        size=(2, 1)), sg.InputText(size=(10, 1))],
        [sg.Text('y0',        size=(2, 1)), sg.InputText(size=(10, 1))],
        [sg.Text('X',         size=(2, 1)), sg.InputText(size=(10, 1))],
        [sg.Text('N',         size=(2, 1)), sg.InputText(size=(10, 1)), sg.Text('N_start',   size=(5, 1)), sg.InputText(size=(10, 1)),\
        sg.Text('N_finish',   size=(6, 1)), sg.InputText(size=(10, 1))],
        [sg.Button('Plot')],
        [sg.Canvas(key='controls_cv')],
        [sg.Text('Figure:')],
        [sg.Column(
            layout=[
                [sg.Canvas(key='fig_cv',
                           # size of the gui
                           size=(800, 300)
                           )]
            ],
            background_color='#DAE0E6',
            pad=(0, 0)
        )],
        [sg.Button('Exit')]
    ]

    window = sg.Window('DE approximation with controls', layout)  # '#e7f4f4'

    while True:
        event, values = window.read()

        if event in (sg.WIN_CLOSED, 'Exit'):
            break

        is_error = 0

        task = values[0]
        is_euler = values[1]
        is_improved_euler = values[2]
        is_runge_kutta = values[3]
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

        X = str(values[6])

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

        if N <= 0:
            if not is_error:
                sg.popup('N should be greater than zero and integer')
                is_error = 1

        if x0 >= X:
            if not is_error:
                sg.popup('x0 must be LESS than X')
                is_error = 1

        if is_error:
            continue

        h = X - x0

        try:
            h = float((X - x0) / N)
        except:
            h = X - x0

        if h > 1:
            is_error = 1
            sg.popup(f'Step must not be greater than 1, current step = {h}')

        if is_error:
            continue

        euler = Euler(h, x0, y0, X, N)
        improved_euler = ImprovedEuler(h, x0, y0, X, N)
        runge_kutta = RungeKutta(h, x0, y0, X, N)

        euler_values = euler.work_method()
        improved_euler_values = improved_euler.work_method()
        runge_kutta_values = runge_kutta.work_method()

        #    return (x_values, runge_kutta_values, exact_values, local_error_value, global_error_values)

        if event == 'Plot':
            check_local_error = 0

            if task == 'GTE for N':
                if not is_euler and not is_improved_euler and not is_runge_kutta:
                    sg.popup('Choose method for which you want to plot GTE')
                    check_local_error = 1

            if task not in ['Approximation', 'LTE', 'GTE for N'] and not check_local_error:
                sg.popup('Unknown task, please change it to the ones that are available')
                check_local_error = 1

            N_start = -1
            N_finish = -1
            if task in ['GTE for N'] and not check_local_error:
                try:
                    N_start = int(values[8])
                except:
                    if not is_error:
                        sg.popup('N_start should be greater than zero and integer')
                        check_local_error = 1

                try:
                    N_finish = int(values[9])
                except:
                    if not is_error:
                        sg.popup('N_finish should be greater than zero and integer')
                        check_local_error = 1

                if N_finish > 1000:
                    if not is_error:
                        sg.popup('N_finish should not be greater than 1000')
                        check_local_error = 1

            if task in ['GTE for N'] and (N_start <= 0 or N_finish <= 0 or N_start > N_finish) and not check_local_error:
                sg.popup('N_start and N_finish should be > 0, and N_start < N_finish')
                check_local_error = 1

            if check_local_error:
                continue

            plt.clf()
            plt.figure(figsize=(10,10))
            fig = plt.gcf()
            DPI = fig.get_dpi()


            title = task
            x_label = 'X'
            y_label = 'Y'

            if task == 'Approximation':
                title = f"Approximation for y' = 3 * y ^ (2 / 3)"

                x_exact = np.linspace(x0, X, 1000)
                y_exact = np.power(np.subtract(x_exact, 1), 3)

                plt.plot(x_exact, y_exact, 'g', label='exact solution')

                if is_euler:
                    plt.plot(euler_values[0], euler_values[1], 'b', label='Euler')

                if is_improved_euler:
                    plt.plot(improved_euler_values[0], improved_euler_values[1], 'y', label='Improved Euler')

                if is_runge_kutta:
                    plt.plot(runge_kutta_values[0], runge_kutta_values[1], 'r', label='Runge-Kutta')

            elif task == 'LTE':
                if is_euler:
                    plt.plot(euler_values[0], euler_values[3], 'b', label='Euler LTE')

                if is_improved_euler:
                    plt.plot(improved_euler_values[0], improved_euler_values[3], 'y', label='Improved Euler LTE')

                if is_runge_kutta:
                    plt.plot(runge_kutta_values[0], runge_kutta_values[3], 'r', label='Runge-Kutta LTE')

            elif task == 'GTE for N':
                # sg.popup("Make sure entered value of 'N' is the end of the range for GTE")
                title = f'GTE for N = [{N_start}, {N_finish}]'

                euler_gte = []
                improved_euler_gte = []
                runge_kutta_gte = []

                n_values = list(range(N_start, N_finish + 1))
                for i in n_values:
                    h_2 = 1
                    try:
                        h_2 = (X - x0) / i
                    except ZeroDivisionError:
                        h_2 = 1

                    euler_2 = Euler(h_2, x0, y0, X, i)
                    improved_euler_2 = ImprovedEuler(h_2, x0, y0, X, i)
                    runge_kutta_2 = RungeKutta(h_2, x0, y0, X, i)

                    euler_values_2 = euler_2.work_method()
                    improved_euler_values_2 = improved_euler_2.work_method()
                    runge_kutta_values_2 = runge_kutta_2.work_method()

                    euler_gte.append(euler_values_2[4])
                    improved_euler_gte.append(improved_euler_values_2[4])
                    runge_kutta_gte.append(runge_kutta_values_2[4])

                if is_euler:
                    plt.plot(n_values, euler_gte, 'b', label='Euler GTE')

                if is_improved_euler:
                    plt.plot(n_values, improved_euler_gte, 'y', label='Improved Euler GTE')

                if is_runge_kutta:
                    plt.plot(n_values, runge_kutta_gte, 'r', label='Runge-Kutta GTE')

                x_label = 'N'
                y_label = 'Max GTE'

            plt.title(f'{title}')
            plt.xlabel(x_label)
            plt.ylabel(y_label)
            plt.grid()
            plt.legend(fontsize=15, loc=0)

            draw_figure_w_toolbar(window['fig_cv'].TKCanvas, fig, window['controls_cv'].TKCanvas)

        else:
            sg.Popup('Ok clicked', keep_on_top=True)

    window.close()


if __name__ == '__main__':
    orig_stdout = sys.stdout
    f = open('output.txt', 'w')
    sys.stdout = f

    gui()

    sys.stdout = orig_stdout
    f.close()
