# -*- coding: utf-8 -*-
import base64
from time import sleep
import pandas as pd
import requests

URL = 'http://www.letpub.com.cn/nsfcfund_search.php?mode=bysubject_s&datakind=excel&currentpage=1'
COOKIES = {'PHPSESSID': 'je9d3no1uvshgkhgr7pcbf3vh3'}

def build_body_data():
    body_list = []
    df = pd.read_excel('/data/dk/dump/中国医学学校.xlsx',engine='openpyxl')
    for college in df['学校名称'].values:
        body_list.append(
            {
                'page': '',
                'organization_main': college,
                'startTime': '2015',
                'endTime': '2023',
                'searchsubmit': True
            }
        )
    return body_list


def write_to_excel(data, file_name):
    print(data)
    response = requests.post(URL, cookies=COOKIES, data=data)

    try:
        result = response.json()['data']
        file_path = f'/data/dk/dump/{file_name}.xls'
        with open(file_path, 'wb') as file:
            result = result.replace('data:application/vnd.ms-excel;base64,', '')
            file.write(base64.b64decode(result))
    except Exception as e:
        print(e)
        print(response.text)


def click():
    url = 'http://www.letpub.com.cn/content/index.php?action=onlinecheck'
    requests.get(url, cookies=COOKIES)


if __name__ == '__main__':

    form_data_list = build_body_data()
    for form_data in form_data_list:
        name = form_data['organization_main'] 
        click()
        write_to_excel(form_data, name)
        print(f'{form_data["organization_main"]} done')
        sleep(10)
