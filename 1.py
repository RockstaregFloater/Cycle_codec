from textwrap import wrap
import numpy as np
from tkinter import *
from random import *
from tkinter import messagebox

def syndroms_and_errors(n, g, Gx, v):
    e, s, c = [], [], []

    for i in range(n + 1):
        e.append([0] * n)
        s.append([0] * n)

    g_x = np.poly1d(list(reversed(g)))
    t = n - 1

    for i in range(1, n + 1):
        e[i][t] = 1
        t -= 1

    for i in range(n + 1):
        e_x = np.poly1d(list((e[i])))
        not_need, s[i] = np.polydiv(e_x, g_x)
        s[i] = list(reversed(s[i].coeffs))
        for j in range(len(s[i])):
            s[i][j] = int(s[i][j])
            if s[i][j] < 0:
                s[i][j] = s[i][j] * -1
            if s[i][j] % 2 == 0 and s[i][j] != 0:
                s[i][j] = 0
            elif s[i][j] % 2 == 1 and s[i][j] != 0:
                s[i][j] = 1

        if len(s[i]) < int(Gx[-1][2]) and i != n:
            for u in range(int(Gx[-1][2]) - len(s[i])):
                s[i].append(0)
        elif i == n:
            for u in range(int(Gx[-1][2]) - len(s[i])):
                s[i].insert(0, 0)

    v_x = np.poly1d(list(reversed(v)))
    not_n, s_x = np.polydiv(v_x, g_x)
    s_x = list(reversed(s_x.coeffs))
    for i in range(len(s_x)):
        s_x[i] = int(s_x[i])
        if s_x[i] < 0:
            s_x[i] = s_x[i] * -1
        if s_x[i] % 2 == 0 and s_x[i] != 0:
            s_x[i] = 0
        elif s_x[i] % 2 == 1 and s_x[i] != 0:
            s_x[i] = 1
        if len(s_x) < int(Gx[-1][2]):
            for u in range(int(Gx[-1][2]) - len(s_x)):
                s_x.append(0)

    for i in range(len(s)):
        if s_x == s[i]:
            e_X = list(reversed(e[i]))
            for j in range(len(e_X)):
                c_x_i = v[j] - e_X[j]
                c.append(c_x_i)

    c_x = np.poly1d(list(reversed(c)))
    i_x, not_n_n = np.polydiv(c_x, g_x)
    i_x = list(reversed(i_x.coeffs))
    for i in range(len(i_x)):
        i_x[i] = int(i_x[i])
        if i_x[i] < 0:
            i_x[i] = i_x[i] * -1
        if i_x[i] % 2 == 0 and i_x[i] != 0:
            i_x[i] = 0
        elif i_x[i] % 2 == 1 and i_x[i] != 0:
            i_x[i] = 1
    if len(i_x) < (int(Gx[-1][2]) + 1):
        for u in range((int(Gx[-1][2]) + 1) - len(i_x)):
            i_x.append(0)

    return i_x

def run():
    window_c = Tk()
    window_c.title("Процесс кодирования")
    window_c.geometry("400x400")

    text = str(text_entry.get())
    if text == "":
        messagebox.showerror("Ошибка", "Пожалуйста, введите текст.")
        window_c.destroy()
        return
    if var.get() == 'polynomial':
        k = 0
        g = g_entry.get()
        Gx = g.split("+")
        n = int(n_entry.get())
        k = n - int(Gx[-1][2])
        if n == '':
            messagebox.showerror("Ошибка", "Пожалуйста, введите текст.")
            window_c.destroy()
            return
        elif n <= k:
            messagebox.showerror("Ошибка", "Пожалуйста, введите n > k.")
            window_c.destroy()
            return
    elif var.get() == 'matrix':
        g = ''
        k = 0
        Gx = []
        g_m = str(m_entry.get())
        if g_m == "":
            messagebox.showerror("Ошибка", "Пожалуйста, введите первую строку матрицы.")
            window_c.destroy()
            return
        n = int(len(g_m))
        for i in range(len(g_m)-1, 0, -1):
            if g_m[i] == '1':
                k = n - i
                break
        for j in range(len(g_m)):
            if g_m[j] == '1':
                Gx.append(f'x^{j}')
        for i in range(len(Gx)):
            g += ''.join(Gx[i]) + '+'
        g = g[:-1]

    text_to_binary = "".join(format(x, "08b") for x in bytearray(text, "utf-8"))
    razbiv = wrap(text_to_binary, k)
    print(text_to_binary)
    if len(razbiv[-1]) < k:
        razbiv[-1] = "0" * (k - len(razbiv[-1])) + razbiv[-1]
    print(razbiv)

    text_info_b = Label(window_c, text='Двоичное представление')
    text_info_b.grid(row=0, column=1)
    out_window_b = Text(window_c, width=40, height=12)
    out_window_b.grid(row=1, column=1)
    for i in range(len(razbiv)):
        out_window_b.insert(END, razbiv[i])
        out_window_b.insert(END, "\n")
    out_window_b.config(state='disabled')

    G = [[0] * n for i in range(k)]
    g_v = G[0]

    for i in range(n):
        for j in range(len(Gx)):
            if i == int(Gx[j][2]):
                G[0][i] = 1
    for i in range(1, k):
        for j in range(n):
            if G[i - 1][0] == 1:
                G[i][1] = 1
            if G[i - 1][j - 1] == 1 and j > 0:
                G[i][j] = 1

    print(G)

    text_info_a = Label(window_c, text='Таблица G')
    text_info_a.grid(row=0, column=0)
    out_window_a = Text(window_c, width=40, height=12)
    out_window_a.grid(row=1, column=0)
    for i in range(len(G)):
        out_window_a.insert(END, G[i])
        out_window_a.insert(END, "\n")
    out_window_a.config(state='disabled')

    a = np.array(G)
    c = razbiv.copy()
    for i in range(len(razbiv)):
        b = [int(x) for x in razbiv[i]]
        c[i] = list(np.matmul(b, a))
    for i in range(len(c)):
        for j in range(n):
            if c[i][j] != 0 and c[i][j] % 2 == 0:
                c[i][j] = 0
            elif c[i][j] != 0 and c[i][j] % 2 != 0:
                c[i][j] = 1

    print(c)
    text_info_c = Label(window_c, text='Закодированные блоки')
    text_info_c.grid(row=0, column=2)
    out_window_c = Text(window_c, width=40, height=12)
    out_window_c.grid(row=1, column=2)
    for i in range(len(c)):
        out_window_c.insert(END, c[i])
        out_window_c.insert(END, "\n")
    out_window_c.config(state='disabled')

    d = c.copy()
    for i in range(len(c)):
        r = 0
        for j in range(len(c[i])):
            if c[i][j] == 1:
                if r > 0:
                    d[i] += f'x^{j}+'
                    r += 1
                else:
                    d[i] = '' + f'x^{j}+'
                    r += 1
        d[i] = d[i][:-1]

    print(d)

    text_info_d = Label(window_c, text='Представление в виде полиномов')
    text_info_d.grid(row=2, column=0)
    out_window_d = Text(window_c, width=40, height=12)
    out_window_d.grid(row=3, column=0)
    for i in range(len(d)):
        out_window_d.insert(END, d[i])
        out_window_d.insert(END, "\n")
    out_window_d.config(state='disabled')

    f = []
    for i in range(len(c)):
        r = c[i]
        f.append(r)
        error_position = randint(0, n - 1)
        for t in range(len(f[i])):
            if t == int(error_position):
                if f[i][t] == 1:
                    f[i][t] = 0
                else:
                    f[i][t] = 1

    print(f)

    text_info_t = Label(window_c, text='Ошибки в блоках')
    text_info_t.grid(row=4, column=0)
    out_window_t = Text(window_c, width=40, height=12)
    out_window_t.grid(row=5, column=0)
    for i in range(len(f)):
        out_window_t.insert(END, f[i])
        out_window_t.insert(END, "\n")
    out_window_t.config(state='disabled')

    e = []
    for i in range(len(f)):
        v = f[i]
        e.append(syndroms_and_errors(n, g_v, Gx, v))

    print(e)

    text_info_e = Label(window_c, text='Декодированные блоки')
    text_info_e.grid(row=2, column=1)
    out_window_e = Text(window_c, width=40, height=12)
    out_window_e.grid(row=3, column=1)
    for i in range(len(e)):
        out_window_e.insert(END, e[i])
        out_window_e.insert(END, "\n")
    out_window_e.config(state='disabled')

    d_str = ''
    for t in range(len(e)):
        for y in range(len(e[t])):
            d_str += ''.join(str(e[t][y]))
    d_str = int(d_str, 2).to_bytes((len(d_str) + 7) // 8, byteorder='big').decode('utf-8')
    print(d_str)

    text_info_f = Label(window_c, text='Декодированная строка')
    text_info_f.grid(row=2,column=2)
    out_window_f = Text(window_c, width=40, height=12)
    out_window_f.grid(row=3,column=2)
    out_window_f.insert(END, d_str)
    out_window_f.config(state='disabled')

    text_info_y = Label(window_c, text='Полином g(x)')
    text_info_y.grid(row=4, column=1)
    out_window_y = Text(window_c, width=40, height=12)
    out_window_y.grid(row=5, column=1)
    out_window_y.insert(END, g)
    out_window_y.config(state='disabled')

def show_entry_fields():
    if var.get() == 'polynomial':
        g_label.pack()
        g_entry.pack()
        n_label.pack()
        n_entry.pack()
        m_label.pack_forget()
        m_entry.pack_forget()
    elif var.get() == 'matrix':
        m_label.pack()
        m_entry.pack()
        n_label.pack_forget()
        n_entry.pack_forget()
        g_label.pack_forget()
        g_entry.pack_forget()

window = Tk()
window.title("Циклические коды")
window.geometry("400x400")

text_label = Label(window, text="Введите текст:")
text_label.pack()

text_entry = Entry(window)
text_entry.pack()

var = StringVar(value='polynomial')

R1 = Radiobutton(window, text="Полином", variable=var, value='polynomial', command=show_entry_fields)
R1.pack()

R2 = Radiobutton(window, text="Матрица", variable=var, value='matrix', command=show_entry_fields)
R2.pack()

g_label = Label(window, text="Образец ввода полинома: x^0+x^1+x^3\nВведите порождающий полином:")

g_entry = Entry(window)

m_label = Label(window, text="Введите первую строку матрицы G, начиная слева:")

m_entry = Entry(window)

n_label = Label(window, text="Введите n:")

n_entry = Entry(window)

run_button = Button(window, text="запуск", command=run)
run_button.pack()

show_entry_fields()

window.mainloop()
