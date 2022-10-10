#!/usr/bin/env python
import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import argparse
from matplotlib.backends.backend_pdf import PdfPages

# #Проверяемый файл
# srcfile
# #Сводная таблица с данными из репортов
# srctable

def Main(srcfile, srctable, Directory=os.path.dirname(os.path.abspath(__file__))):
    #Основной dataframe задается из таблицы
    dfall = tableread(srctable)

    #Выполняется подгрузка файла в таблицу
    if srcfile in dfall.index:
      # Искомый файл есть в базе, вносить новый репорт в базу не нужно
      pass
    else:
      # Искомый файл в базе отсутствует, для проверки необходимо внести его в базу
      dfall = pd.concat([dfall, get_file_df(srcfile)]).drop_duplicates()
      dfall.to_csv('Maintable2.tsv', sep='\t')

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
    filename = "multi"
    #Сохранение графиков в виде pdf файла
    save_multi_png(filename)
    plt.interactive(False)
    #Вызов функции проверки файла
    checkSF(dfall, srcfile, Directory)

def tableread(table):
    """
    Читает .csv/.tsv файл и формирует общий фрейм
    :param table(str): название таблицы
    :return df(pandas.DataFrame): общий фрейм из сводной таблицы
    """
    df = pd.read_csv(table, sep='\t')
    df.set_index('Unnamed: 0', inplace=True)
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
def get_file_df(filename):
  """
  На выделяет данны метрик из файла отчета в 
  pandas DataFrame
  :param filename(str): имя файла
  :return(pandas.DataFrame): фрейм отдельного репорта
  """
  lstdf = DFSplit(aquire(filename))
  df1 = lstdf[0].to_frame().T[['Clean read rate', 'Mean Coverage', 'Raw reads']]
  df2 = lstdf[1].to_frame().T[['Total filtered', 'Cutted adapter']]
  depdis10 = lstdf[3].to_frame().T.add_suffix('_chr10')
  depdisX = lstdf[4].to_frame().T.add_suffix('_chrX')
  depdisY = lstdf[5].to_frame().T.add_suffix('_chrY')
  df7 = lstdf[6].to_frame().T
  df8 = lstdf[7].to_frame().T.add_suffix('_chr10')
  df9 = lstdf[8].to_frame().T.add_suffix('_chrX')
  df10 = lstdf[9].to_frame().T.add_suffix('_chrY')
  lstmdf = [df1, df2, depdis10, depdisX, depdisY, df7, df8, df9, df10]
  for df in lstmdf:
    df.index = [filename]
  outdf = pd.concat(lstmdf, axis=1).replace({'\%': ''}, regex=True).astype(float).dropna(how='any')
  return outdf
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

def save_multi_png(filename):
    fig_nums = plt.get_fignums()
    figs = [plt.figure(n) for n in fig_nums]
    for i, fig in enumerate(figs):
        fig.savefig("{0}_{1}.png".format(i, filename))
def checkSF(maindf, checkstr, Directory):
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
    if data[['dbSNP rate_chr10']].iat[0, 0] < 99:
        q1=' < 99 '
        q2=' '
        q3=' '
        out = "Fail"
    elif data[['Mean Coverage']].iat[0, 0] < 15:
        q1=' '
        q2=' < 15 '
        q3=' '
        out = "Fail"
    elif loval>data[['dbSNP TITV_chr10']].iat[0, 0]>hival:
        q1=' '
        q2=' '
        q3=' outside boxplot whiskers '
        out = "Fail"

    else:
        q1 = ' > 99 '
        q2 = ' > 15 '
        q3 = ' within boxplot whiskers '
        out = "Success"


    # print(data[['Mean Coverage']].iat[0, 0] + ' < 15 ')
    # print(data[['dbSNP rate_chr10']].iat[0, 0] + ' в пределах "усов"')
    d = {'dbSNP rate_chr10':str(data[['dbSNP rate_chr10']].iat[0, 0])+q1, 'Mean Coverage':str(data[['Mean Coverage']].iat[0, 0])+q2,
         'dbSNP TITV_chr10':str(data[['dbSNP TITV_chr10']].iat[0, 0])+q3, 'Name':checkstr, 'Result':out}
    finalscore = pd.Series(data=d, name='Summary')
    finalscore.to_csv(Directory+'finalscore.txt', sep='\t')
    return finalscore

parser = argparse.ArgumentParser()
parser.add_argument("srcf", help="input file name to check")
parser.add_argument("srct", help="input unified table location")
parser.add_argument("Dir", help="output directory")
args = parser.parse_args()

if __name__ == '__main__':
    Main(srcfile=args.srcf, srctable=args.srct,
         Directory=args.Dir)
    # plt.show()
