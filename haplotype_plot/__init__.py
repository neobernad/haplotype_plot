# -*- coding: utf-8 -*-
import logging

logging.basicConfig(format='%(asctime)s %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s',
                    datefmt='%Y-%m-%d:%H:%M:%S')
version_info = (1, 1, 1)
version = '.'.join(str(c) for c in version_info)
# https://codereview.stackexchange.com/questions/129683/storing-a-version-number-in-a-python-project
