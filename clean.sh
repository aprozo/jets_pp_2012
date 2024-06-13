echo "cleaning that directory"
rm *session.xml
rm -r run12data.package
rm *.zip
rm schedTemplateExp.xml
rm -r sched*.package
rm *.dataset

echo "additional out/ report/ csh/"
rm submit/scheduler/csh/*
rm submit/scheduler/list/*
rm submit/scheduler/report/*
rm submit/log/*
