
import cgi
import os
import sys
sys_path = os.path.dirname(os.getcwd())
src_path = sys_path + "/infection_risk_calculator/src"
sys.path.insert(1, src_path)
from calculator import * 
src_path = sys_path + "/infection_risk_calculator"
sys.path.insert(1, src_path)
form = cgi.FieldStorage()

n_occ = form.getvalue("nocc")
t = form.getvalue("t")
r_id = str(form.getvalue("rid"))
act_input = str(form.getvalue('act'))

activity = 'nan'
exp_activity = 'nan'
if act_input == 'Lecture':
    activity = "standing"
    exp_activity = "speaking"
if activity == 'Studying':
    activity = 'standing'
    exp_activity = 'whispering'
if activity == 'Singing':
    activity = 'standing'
    exp_activity = 'singing'
if activity == 'Social Event':
    activity = 'light_exercise'
    exp_activity = 'speaking'

ir = infection_risk(t, r_id, n_occ, activity, exp_activity, sys_path + '/data/raw/rm.csv')

print('<html>')
print('<body>')
print('<h1> Infection Risk: %s <h1>'%ir)
print('<body>')
print('<html>')