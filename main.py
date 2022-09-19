import re
import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

# #Проверяемый файл
# srcfile
# #Сводная таблица с данными из репортов
# srctable
# #Директория с репортами
# Directory
# #Режимы

#readfolder - программа читает файлы репортов из заданной папки
#readtable - программа чиатет сводную таблицу
Mode = 'readtable'

def Main(srcfile, srctable, Mode, Directory=os.path.dirname(os.path.abspath(__file__))):
    #Программа проверяет режим
    if Mode == 'readtable':
        #Основной dataframe задается из таблицы
        dfall = tableread(srctable)

    elif Mode == 'readfolder':
        # Основной dataframe задается из папки
        dfall = fdclean(Directory)
        # Данные сохраняются в виде таблицы
        dfall.to_csv('Maintable2.tsv', sep='\t')
    else:
        raise ValueError('Mode not selected')

    #Задаются датасеты для построения карт коррелляции
    corrm = dfall[dfall['Total SNP_chrX'] < 90000].corr().replace({np.NAN: 0})
    corrf = dfall[dfall['Total SNP_chrX'] > 90000].corr().replace({np.NAN: 0})

    #Определяются параметры фигур
    plt.rcParams["figure.figsize"] = [10.00, 10.00]
    plt.rcParams["figure.autolayout"] = True

    #Отрисовка графиков
    plot_QCres(dfall[['Clean read rate']], dfall[['Mean Coverage']], dfall[['Raw reads']], srcfile)
    plot_DetTrim(dfall[['Total filtered']], dfall[['Cutted adapter']], srcfile)
    plot_DepthDistr(dfall[['Coverage(>=1X)_chr10', 'Coverage(>=5X)_chr10', 'Coverage(>=10X)_chr10', 'Coverage(>=20X)_chr10']],
                    dfall[['Coverage(>=1X)_chrX', 'Coverage(>=5X)_chrX', 'Coverage(>=10X)_chrX', 'Coverage(>=20X)_chrX']],
                    dfall[['Coverage(>=1X)_chrY', 'Coverage(>=5X)_chrY', 'Coverage(>=10X)_chrY', 'Coverage(>=20X)_chrY']])
    plot_AlRes(dfall[['Mapping rate']], dfall[['Mismatch rate']], srcfile)
    plot_VarStat(dfall[['Total SNP_chr10']], dfall[['dbSNP rate_chr10']], dfall[['Total SNP_chrX']],
                 dfall[['dbSNP rate_chrX']], dfall[['Total SNP_chrY']], dfall[['dbSNP rate_chrY']], srcfile)
    plot_CorrFM(corrf, corrm)

    #Имя pdf файла с графиками
    filename = "multi.pdf"
    #Сохранение графиков в виде pdf файла
    save_multi_image(filename)
    #Вызов функции проверки файла
    checkSF(dfall, srcfile)

def tableread(table):
    """
    Читает .csv/.tsv файл и формирует общий фрейм
    :param table(str): название таблицы
    :return df(pandas.DataFrame): общий фрейм из сводной таблицы
    """
    df = pd.read_csv(table, sep='\t')
    df.set_index('Unnamed: 0', inplace=True)
    return df
def aquire(file):
    """
    Читает html файл репорта и формирует фрейм
    на основе таблиц
    :param file(str): имя файла
    :return(pandas.DataFrame): фрейм отдельного репорта
    """
    data = pd.read_html(file)
    out = DFFormat(data)
    return out
def prune(df):
    """
    Удаляет символы %, конвертирует значения
    общего фрейма в float, устраняет бесконечные значения
    :param df(pandas.DataFrame): общий фрейм
    :return: df(pandas.DataFrame): общий фрейм
    """
    df = df.replace({'\%': ''}, regex=True).astype(float).dropna(how='any')
    df.replace([np.inf, -np.inf], np.nan, inplace=True).dropna(how='any')
    return df
def DFFormat(dfs):
    """
    Устанавливает индекс строк
    :param dfs(pandas.DataFrame): фрейм из фреймов
    :return: dfs(pandas.DataFrame): фрейм из фреймов
    """
    for df in dfs:
        df.set_index('Metric', inplace=True)
    return(dfs)
def DFSplit(dfs):
    """
    Разделяет исходные данные из папки (фрейм из фреймов)
    на отдельные фреймы
    :param dfs(pandas.DataFrame): фрейм из фреймов
    :return: dflst(list) список из разделенных фремов
    """
    dflst = []
    for df in dfs:
        cols = list(df)
        for col in cols:
            ndf = df[col]
            dflst.append(ndf)
    return dflst
def TableBuilder(dir):
    """
    Собирает фрейм из фреймов на основе списка фреймов,
    выбирает только файлы подходящие под паттерн названия
    :param dir(str): путь к расположению папки с репортами
    :return: df(pandasDataFrame) фрейм из фреймов
    """
    #Создает список из названий файлов репортов
    lst = os.listdir(dir)
    #ре паттерн выбора правильных названий
    a = re.compile(r'.{11}[.].{5}[.]WGS[.]report[.]html')
    newlist = []
    chnum = 0
    for file in lst:
        if a.match(file) != None:
            nfile = dir+file
            chnum+=1
            lst = DFSplit(aquire(nfile))
            for item in lst:
                item.name = file
            newlist.append(lst)
    df = pd.DataFrame(newlist)
    return df
def fdclean(dir):
    """
    Отбрасывает ненужные метрики. добавляет суффиксы к метрикам с одинаковым названием
    Вытаскивает отдельные фреймы из фрейма фреймов, собирает в общий датафрейм
    :param dir(str): путь к папке с репортами
    :return: bigdf(pandas.DataFrame) общий датафрейм
    """
    df = TableBuilder(dir)
    df1 = pd.DataFrame(item for item in df[0])[['Clean read rate', 'Mean Coverage', 'Raw reads']]
    df2 = pd.DataFrame(item for item in df[1])[['Total filtered', 'Cutted adapter']]
    depdis10 = pd.DataFrame(item for item in df[3]).add_suffix('_chr10')
    depdisX = pd.DataFrame(item for item in df[4]).add_suffix('_chrX')
    depdisY = pd.DataFrame(item for item in df[5]).add_suffix('_chrY')
    df7 = pd.DataFrame(item for item in df[6])
    df8 = pd.DataFrame(item for item in df[7]).add_suffix('_chr10')
    df9 = pd.DataFrame(item for item in df[8]).add_suffix('_chrX')
    df10 = pd.DataFrame(item for item in df[9]).add_suffix('_chrY')
    bigdf = prune(pd.concat([df1, df2, depdis10, depdisX, depdisY, df7, df8, df9, df10], axis=1))
    return bigdf
def plot_QCres(df1, df2, df3, x):
    fig, axs = plt.subplots(1, 3)
    sns.boxplot(ax=axs[0], data=df1, fliersize=0, width=0.3)
    sns.stripplot(ax=axs[0], data=df1, color="darkorange", jitter=0.12, dodge=True, size=1.5)
    plot_Val(df1, axs[0], x)
    sns.boxplot(ax=axs[1], data=df2, fliersize=0, width=0.3)
    sns.stripplot(ax=axs[1], data=df2, color="darkorange", jitter=0.12, dodge=True, size=1.5)
    plot_Val(df2, axs[1], x)
    sns.boxplot(ax=axs[2], data=df3, fliersize=0, width=0.3)
    sns.stripplot(ax=axs[2], data=df3, color="darkorange", jitter=0.12, dodge=True, size=1.5)
    plot_Val(df3, axs[2], x)
    return fig
def plot_DetTrim(df1, df2, x):
    fig, axs = plt.subplots(1, 2)
    sns.boxplot(ax=axs[0], data=df1, fliersize=0, width=0.3)
    sns.stripplot(ax=axs[0], data=df1, color="darkorange", jitter=0.12, dodge=True, size=1.5)
    plot_Val(df1, axs[0], x)
    sns.boxplot(ax=axs[1], data=df2, fliersize=0, width=0.3)
    sns.stripplot(ax=axs[1], data=df2, color="darkorange", jitter=0.12, dodge=True, size=1.5)
    plot_Val(df2, axs[1], x)
    return fig
def plot_DepthDistr(df1, df2, df3):
    fig, axs = plt.subplots(1, 3)
    sns.boxplot(ax=axs[0], data=df1, fliersize=0, width=0.3)
    sns.stripplot(ax=axs[0], data=df1, color="darkorange", jitter=0.12, dodge=True, size=1.5)
    sns.boxplot(ax=axs[1], data=df2, fliersize=0, width=0.3)
    sns.stripplot(ax=axs[1], data=df2, color="darkorange", jitter=0.12, dodge=True, size=1.5)
    sns.boxplot(ax=axs[2], data=df3, fliersize=0, width=0.3)
    sns.stripplot(ax=axs[2], data=df3, color="darkorange", jitter=0.12, dodge=True, size=1.5)
    for ax in axs:
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    return fig
def plot_AlRes(df1, df2, x):
    fig, axs = plt.subplots(1, 2)
    sns.boxplot(ax=axs[0], data=df1, fliersize=0, width=0.3)
    sns.stripplot(ax=axs[0], data=df1, color="darkorange", jitter=0.12, dodge=True, size=1.5)
    plot_Val(df1, axs[0], x)
    sns.boxplot(ax=axs[1], data=df2, fliersize=0, width=0.3)
    sns.stripplot(ax=axs[1], data=df2, color="darkorange", jitter=0.12, dodge=True, size=1.5)
    plot_Val(df2, axs[1], x)
    return fig
def plot_VarStat(df1, df2, df3, df4, df5, df6, x):
    fig, axs = plt.subplots(3, 2)
    sns.boxplot(ax=axs[0][0], data=df1, fliersize=0, width=0.3)
    sns.stripplot(ax=axs[0][0], data=df1, color="darkorange", jitter=0.12, dodge=True, size=1.5)
    plot_Val(df1, axs[0][0], x)
    sns.boxplot(ax=axs[0][1], data=df2, fliersize=0, width=0.3)
    sns.stripplot(ax=axs[0][1], data=df2, color="darkorange", jitter=0.12, dodge=True, size=1.5)
    plot_Val(df2, axs[0][1], x)
    sns.boxplot(ax=axs[1][0], data=df3, fliersize=0, width=0.3)
    sns.stripplot(ax=axs[1][0], data=df3, color="darkorange", jitter=0.12, dodge=True, size=1.5)
    plot_Val(df3, axs[1][0], x)
    sns.boxplot(ax=axs[1][1], data=df4, fliersize=0, width=0.3)
    sns.stripplot(ax=axs[1][1], data=df4, color="darkorange", jitter=0.12, dodge=True, size=1.5)
    plot_Val(df4, axs[1][1], x)
    sns.boxplot(ax=axs[2][0], data=df5, fliersize=0, width=0.3)
    sns.stripplot(ax=axs[2][0], data=df5, color="darkorange", jitter=0.12, dodge=True, size=1.5)
    plot_Val(df5, axs[2][0], x)
    sns.boxplot(ax=axs[2][1], data=df6, fliersize=0, width=0.3)
    sns.stripplot(ax=axs[2][1], data=df6, color="darkorange", jitter=0.12, dodge=True, size=1.5)
    plot_Val(df6, axs[2][1], x)
    return fig
def plot_CorrFM(df1, df2):
    sns.clustermap(data=df1, figsize=(10, 10), cmap='magma').fig.suptitle('Female correlation matrix')
    sns.clustermap(data=df2, figsize=(10, 10), cmap='magma').fig.suptitle('Male correlation matrix')
def plot_Val(df, ax, value):
    mask = df.index.str.startswith(value, na=False)
    col = df.columns
    sns.stripplot(data=df[mask][col], ax = ax,
               marker='x', s=15, linewidth=2, color='red', jitter=0)
def save_multi_image(filename):
    pp = PdfPages(filename)
    fig_nums = plt.get_fignums()
    figs = [plt.figure(n) for n in fig_nums]
    for fig in figs:
        fig.savefig(pp, format='pdf')
    pp.close()
def checkSF(maindf, checkstr):
    """
    Проверяет заданный репорт на вхождение
    в "усы" боксплота по метрике dbSNP TITV_chr10
    и по пороговым значением метрик dbSNP rate_chr10 и Mean Coverage
    :param maindf(pandas.DataFrame): общий датафрейм
    :param checkstr(str): название проверяемого файла
    :return: finalscore(pandas.DataFrame) датафрейм из рассматриваемых характеристик и финальной оценки
    """
    mask = maindf.index.str.startswith(checkstr, na=False)
    data = maindf[mask][['dbSNP rate_chr10', 'Mean Coverage', 'dbSNP TITV_chr10']]
    Q1 = maindf['dbSNP TITV_chr10'].quantile(0.25)
    Q3 = maindf['dbSNP TITV_chr10'].quantile(0.75)
    IQR = Q1-Q3
    loval = Q1 - 1.5 * IQR
    hival = Q3 + 1.5 * IQR
    if data[['dbSNP rate_chr10']].iat[0, 0] < 99 \
            or data[['Mean Coverage']].iat[0, 0] < 15 \
            or loval>data[['dbSNP TITV_chr10']].iat[0, 0]>hival:
        out = "Fail"
    else:
        out = "Success"
    d = {'dbSNP rate_chr10':data[['dbSNP rate_chr10']].iat[0, 0], 'Mean Coverage':data[['Mean Coverage']].iat[0, 0],
         'dbSNP TITV_chr10':data[['dbSNP TITV_chr10']].iat[0, 0], 'Name':checkstr, 'Result':out}
    finalscore = pd.Series(data=d, name='')
    finalscore.to_csv('finalscore.txt', sep='\t')
    return finalscore


if __name__ == '__main__':
    plt.interactive(False)
    Main(srcfile='DZM00000028.16382.WGS.report.html', srctable='Maintable2.tsv', Mode='readtable',
         Directory='C:/Users/Valery.LAPTOP-FAILHIKM/Desktop/reports/')
    # plt.show()