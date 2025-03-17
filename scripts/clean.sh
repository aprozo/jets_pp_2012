echo "cleaning that directory"
rm *session.xml
rm -r *.package/

rm *.zip
rm schedTemplateExp.xml
rm sched*.package
rm *.dataset

echo "additional out/ report/ csh/"
rm -rf submit/scheduler/csh/
rm -rf submit/scheduler/list/
rm -rf submit/scheduler/report/
rm -rf submit/scheduler/gen
rm -rf submit/log/

mkdir -p submit/scheduler/gen
mkdir -p submit/scheduler/csh
mkdir -p submit/scheduler/list
mkdir -p submit/scheduler/report
mkdir -p submit/log/
