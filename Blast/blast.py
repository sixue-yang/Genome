# _author_ = 'sixueyang'
# _date_ = 2023/5/23 9:58
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import multiprocessing as mp
import threading
import queue
import sys
import os
import gzip
import sqlite3
from collections import defaultdict
import pandas as pd
import numpy as np
import yaml
import argparse


# blast线程控制
def blast_worker(model,input_queue, output_queue,thread):
    # 启动子进程，将输入序列从input_queue获取，执行比对，将结果放入output_queue中
    while True:
        try:
            # 设置队列等待时间，若超过1s钟则退出
            record_id,record_seq = input_queue.get(timeout=1)
        except:
            break
        if model == 'NT':
            result = blastn_sequence(str(record_seq),thread)
        else:
            result = blastp_sequence(str(record_seq),thread)
        for i in result[0].split('\n'):
            if i:
                tag = i.strip().split('\t')
                tag[0] = record_id
                output_queue.put(tag)

# blast 核酸比对模块
def blastn_sequence(record,num_threads):
    # 比对单个序列
    blastn_cline = NcbiblastpCommandline(cmd=blastn, query='-', db=NT_blast_db, outfmt=6, evalue=0.00001,num_threads=num_threads)
    result = blastn_cline(stdin=str(record))
    return result

# blast 蛋白比对模块
def blastp_sequence(record,num_threads):
    # 比对单个序列
    blastn_cline = NcbiblastnCommandline(cmd=blastp, query='-', db=NR_blast_db, outfmt=6, evalue=0.00001,num_threads=num_threads)
    result = blastn_cline(stdin=str(record))
    return result


# sqlite 数据库匹配
def search_contig_values(model,db_path,contig_ids):
    # 连接到SQLite3数据库
    conn = sqlite3.connect(db_path)

    # 创建游标对象
    cursor = conn.cursor()

    # 构建SQL语句
    placeholders = ','.join(['?' for _ in range(len(contig_ids))])
    sql = f"SELECT Contig_id,Species,Genus,Kingdom FROM {model} WHERE Contig_id IN ({placeholders})"

    # 执行SQL语句并获取结果
    results = cursor.execute(sql, contig_ids).fetchall()

    # 关闭游标和数据库连接
    cursor.close()
    conn.close()

    return results


# 转化搜索结果数据结构
def trun_match_result(results):
    species_dic = {}
    genus_dic = {}
    kingdom_dic = {}
    for Contig_id,Species,Genus,Kingdom in results:
        species_dic[Contig_id] = Species
        genus_dic[Contig_id] = Genus
        kingdom_dic[Contig_id] = Kingdom
    return species_dic,genus_dic,kingdom_dic





# 设置参数
def get_args():
    # 设置命令行参数解析器
    parser = argparse.ArgumentParser(description="NT、NR库 blast比对")
    parser.add_argument("-i", "--input", help="输入文件", required=True)
    parser.add_argument("-dt", "--datatype", type=str, help="数据类型 fasta 或 fastq", choices=['fasta','fastq'],default='fasta')
    parser.add_argument("-m", "--model", type=str, help="比对的数据库类型 NT 或 NR",choices=['NT','NR'],default='NT')
    parser.add_argument("-t", "--thread", type=int, help="比对使用线程数", default=12)
    parser.add_argument("-o", "--outfile", type=str, help="输出文件",required=True)

    return parser.parse_args()





def main(model,data_type,thread:int,input_file:str,out_file:str):

    # 读取输入文件
    input_queue = queue.Queue()
    if data_type == 'fasta':
        if input_file.endswith('gz'):
            with gzip.open(input_file,"rt") as in_handle:
                for title, seq in SimpleFastaParser(in_handle):
                    input_queue.put((title, seq))
        else:
            with open(input_file,'r') as in_handle:
                for title, seq in SimpleFastaParser(in_handle):
                    input_queue.put((title, seq))
    elif data_type == 'fastq':
        if input_file.endswith('gz'):
            with gzip.open(input_file,"rt") as in_handle:
                for title, seq, qual in FastqGeneralIterator(in_handle):
                    input_queue.put((title, seq))
        else:
            with open(input_file,'r') as in_handle:
                for title, seq, qual in FastqGeneralIterator(in_handle):
                    input_queue.put((title, seq))

    else:
        sys.exit('数据类型存在问题，请重新输入数据类型 fasta or fastq')
    # 创建多个进程来执行blast比对
    output_queue = queue.Queue()
    workers = [threading.Thread(target=blast_worker, args=(model, input_queue, output_queue,thread)) for _ in
               range(2)]
    # 开启线程任务
    for worker in workers:
        worker.start()

    # 等待所有进程执行完毕
    for worker in workers:
        worker.join()

    # 存储blast比对 m8 结果
    result_dic = defaultdict(list)

    header = ['Query_ID', 'Subject_ID', 'Percent_identity', 'Alignment_length', 'Number_of_mismatches',
              'Number_of_gap_openings', 'Query_start_position', 'Query_end_position', 'Subject_start_position',
              'Subject_end_position', 'E_value', 'Bit_score']

    # 从结果队列中获取比对结果，并输出
    while not output_queue.empty():
        result = output_queue.get()
        for index, key in enumerate(header):
            result_dic[key].append(result[index])

    if model == 'NT':
        db_path = NT_contig_db
    else:
        db_path = NR_contig_db
    # 匹配搜索结果
    result_df = pd.DataFrame(result_dic)
    search_contig_list = list(set(result_df['Subject_ID'].to_list()))
    results = search_contig_values(model, db_path, search_contig_list)
    species_dic, genus_dic, kingdom_dic = trun_match_result(results)

    result_df['Species'] = result_df['Subject_ID'].map(species_dic)
    result_df['Genus'] = result_df['Subject_ID'].map(genus_dic)
    result_df['Kingdom'] = result_df['Subject_ID'].map(kingdom_dic)

    result_df.to_csv(out_file, sep='\t', index=False)


if __name__ == '__main__':

    # 导入参数
    args = get_args()
    input_file = args.input
    datatype = args.datatype
    model = args.model
    thread = args.thread
    outfile = args.outfile

    # 加载配置文件，获取配置信息
    with open(os.path.join(os.path.dirname(sys.argv[0]), 'db.config'), 'r', encoding='utf-8') as config_obj:
        config = yaml.load(config_obj, Loader=yaml.FullLoader)

        # 比对数据库路径
        NT_blast_db = config['NT_blast_db']
        NR_blast_db = config['NR_blast_db']
        # contig对应关系数据库路径
        NT_contig_db = config['NT_contig_db']
        NR_contig_db = config['NR_contig_db']
        # 软件路径
        blastn = config['blastn']
        blastp = config['blastp']

        # 开始进行blast比对主流程
        main(model,datatype,thread,input_file,outfile)









