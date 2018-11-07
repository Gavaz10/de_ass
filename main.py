import sys
import math
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QHBoxLayout, QSizePolicy, QPushButton, QLineEdit, \
    QLabel, QComboBox
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

e = 2.718281

x_graph = []
y_graph = []
y_graph_error = []
x_appr = []
y_appr = []


def runge_kutta(f, y0, x0, X, h, ap=False):
    x, y = x0, y0
    error_max = -math.inf
    while x <= X:
        k1 = f(x, y)
        k2 = f(x + 0.5 * h, y + 0.5 * h * k1)
        k3 = f(x + 0.5 * h, y + 0.5 * h * k2)
        k4 = f(x + h, y + h * k3)
        y += h * (k1/6 + k2/3 + k3/3 + k4/6)
        if not ap:
            x_graph.append(x)
            y_graph.append(y)
        else:
            error_max = max(error_max, calc_error(y, x0, y0, x))
        x += h
    if ap:
        x_appr.append((X-x0)/h)
        y_appr.append(error_max)


def improved_euler(f, y0, x0, X, h, ap=False):
    t,y = x0,y0
    error_max = -math.inf
    while t <= X:
        k1 = f(t, y)
        k2 = f(t + h, y + h * k1)
        y += h*(k1 + k2)/2
        if not ap:
            x_graph.append(t)
            y_graph.append(y)
        else:
            error_max = max(error_max, calc_error(y, x0, y0, t))
        t += h
    if ap:
        x_appr.append((X-x0)/h)
        y_appr.append(error_max)


def euler(f, y0, x0, X, h, ap=False):
    t,y = x0,y0
    error_max = -math.inf
    while t <= X:
        y += h * f(t,y)
        if not ap:
            x_graph.append(t)
            y_graph.append(y)
        else:
            error_max = max(error_max, calc_error(y, x0, y0, t))
        t += h
    if ap:
        x_appr.append((X-x0)/h)
        y_appr.append(error_max)


def exact_solution(y0, x0, X, h, ap=False):
    t = x0
    error_max = -math.inf
    while t <= X:
        c = y0 / (e ** (-math.sin(x0))) - x0
        y = (t + c) * e ** (-math.sin(t))
        if not ap:
            x_graph.append(t)
            y_graph.append(y)
        else:
            error_max = max(error_max, calc_error(y, x0, y0, t))
        t += h
    if ap:
        x_appr.append((X-x0)/h)
        y_appr.append(error_max)


def calc_error(num_ans, x0, y0, x):
    c = y0 / (e**(-math.sin(x0))) - x0
    y = (x + c) * e**(-math.sin(x))#get exact answer
    return abs(num_ans - y)


def func(x, y):
    return e**(-math.sin(x)) - y*math.cos(x)


class App(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()

    def on_click(self):
        x_graph.clear()
        y_graph.clear()
        y_graph_error.clear()
        x_appr.clear()
        y_appr.clear()
        y0 = float(self.y0.text())
        x0 = float(self.x0.text())
        x = float(self.x.text())
        step = float(self.step.text())
        n0 = int(self.n0.text())
        nf = int(self.nf.text())
        if self.method_box.currentIndex() == 0:
            euler(func, y0, x0, x, step)
        elif self.method_box.currentIndex() == 1:
            improved_euler(func, y0, x0, x, step)
        elif self.method_box.currentIndex() == 2:
            runge_kutta(func, y0, x0, x, step)
        elif self.method_box.currentIndex() == 3:
            exact_solution(y0, x0, x, step)
        for z in range(len(x_graph)):
            i = x_graph[z]
            j = y_graph[z]
            y_graph_error.append(calc_error(j, x0, y0, i))
        for i in range(n0, nf + 1):
            if self.method_box.currentIndex() == 0:
                euler(func, y0, x0, x, (x - x0)/i, ap=True)
            elif self.method_box.currentIndex() == 1:
                improved_euler(func, y0, x0, x, (x - x0)/i, ap=True)
            elif self.method_box.currentIndex() == 2:
                runge_kutta(func, y0, x0, x, (x - x0)/i, ap=True)
            elif self.method_box.currentIndex() == 3:
                exact_solution(y0, x0, x, (x - x0)/i, ap=True)
        self.graph.plot(x_graph, y_graph)
        self.errors.plot(x_graph, y_graph_error)
        self.appr_errors.plot(x_appr, y_appr)


    def initUI(self):
        vbox = QVBoxLayout()

        hbox_graphs = QHBoxLayout()
        hbox_gui = QHBoxLayout()

        self.graph = PlotCanvas(width=5, height=4, text="Graph of solutions")
        self.errors = PlotCanvas(width=5, height=4, text="Graph of errors")
        self.appr_errors = PlotCanvas(width=5, height=4, text="Approximation errors")

        self.method_box = QComboBox()
        self.method_box.addItem("Euler's method")
        self.method_box.addItem("Improved Euler’s method")
        self.method_box.addItem("Runge-Kutta method")
        self.method_box.addItem("Exact solution")
        y0_text = QLabel("y0")
        self.y0 = QLineEdit("1")
        x0_text = QLabel("x0")
        self.x0 = QLineEdit("0")
        x_text = QLabel("X")
        self.x = QLineEdit("9.3")
        step_text = QLabel("Step size")
        self.step = QLineEdit("0.01")
        n0_text = QLabel("n0")
        self.n0= QLineEdit("75")
        nf_text = QLabel("nf")
        self.nf = QLineEdit("400")

        button = QPushButton('Solve')
        button.clicked.connect(self.on_click)

        hbox_graphs.addWidget(self.graph)
        hbox_graphs.addWidget(self.errors)
        hbox_graphs.addWidget(self.appr_errors)

        hbox_gui.addWidget(self.method_box)
        hbox_gui.addWidget(y0_text)
        hbox_gui.addWidget(self.y0)
        hbox_gui.addWidget(x0_text)
        hbox_gui.addWidget(self.x0)
        hbox_gui.addWidget(x_text)
        hbox_gui.addWidget(self.x)
        hbox_gui.addWidget(step_text)
        hbox_gui.addWidget(self.step)
        hbox_gui.addWidget(n0_text)
        hbox_gui.addWidget(self.n0)
        hbox_gui.addWidget(nf_text)
        hbox_gui.addWidget(self.nf)
        hbox_gui.addWidget(button)
        hbox_gui.addStretch()

        vbox.addLayout(hbox_graphs)
        vbox.addLayout(hbox_gui)

        self.setLayout(vbox)
        self.resize(1200, 400)
        self.setWindowTitle('de practicum 11th variant e^(-sin(x)) - y*cos(x)')
        self.show()


class PlotCanvas(FigureCanvas):
    def __init__(self, parent=None, width=5, height=4, dpi=100, text="lul"):
        self.text = text
        self.fig = Figure(figsize=(width, height), dpi=dpi)

        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                                   QSizePolicy.Expanding,
                                   QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    def plot(self, x, y):
        self.fig.clear()
        self.axes = self.fig.add_subplot(111)
        self.axes = self.figure.add_subplot(111)
        self.axes.plot(x, y)
        self.axes.set_title(self.text)
        self.draw()


if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = App()
    sys.exit(app.exec_())