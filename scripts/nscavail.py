#!/usr/bin/env python3

#figure out of an NSC id is actually available by scraping the answer
#interestingly, the NCI people don't seem to be interested in provided a more
#efficient interface for extracting this information

import sys,urllib.request,urllib.error,urllib.parse,urllib.request,urllib.parse,urllib.error,re
import mechanize
from bs4 import BeautifulSoup as bs

def nscavail(id):  #id should be a string
    id = re.sub(r'NSC-?', '', id)
    url = 'https://dtp.cancer.gov/RequestCompounds/forms/order.xhtml'
    browser = mechanize.Browser()
    browser.set_handle_robots(False)
    browser.open(url,timeout=30)
    browser.select_form('orderForm')
    browser.submit('orderForm:j_idt11')
    browser.select_form('orderForm')
    browser['orderForm:nsc'] = id
    browser['orderForm:amt'] = '1'
    browser.submit('orderForm:j_idt30')
    soup = bs(browser.response().read(),"html5lib")
    
    msgs = soup.find_all(class_='ui-messages-info-summary') + soup.find_all(class_='ui-messages-error-summary')
    return len(msgs) == 0

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Need NSC id as argument")
    else:
        id = sys.argv[1]
        print(nscavail(id))
