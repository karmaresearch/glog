/*
 * Copyright 2021 Jacopo Urbani
 *
 * Licensed to the Apache Software Foundation (ASF) under one
 * or more contributor license agreements.  See the NOTICE file
 * distributed with this work for additional information
 * regarding copyright ownership.  The ASF licenses this file
 * to you under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing,
 * software distributed under the License is distributed on an
 * "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
 * KIND, either express or implied.  See the License for the
 * specific language governing permissions and limitations
 * under the License.
 **/


#include <Python.h>
#include <iostream>
#include <vector>

#include <vlog/wizard.h>
#include <python/glog.h>

#include <python/glog.h>

/*** Methods ***/
static PyObject * wizard_new(PyTypeObject *type, PyObject *args, PyObject *kwds);
static int wizard_init(glog_Wizard *self, PyObject *args, PyObject *kwds);
static void wizard_dealloc(glog_Wizard* self);
static PyObject *wizard_rewrite_program(PyObject *self, PyObject *args);


static PyMethodDef Wizard_methods[] = {
    {"rewrite_program", wizard_rewrite_program, METH_VARARGS, "Rewrite a program with magic sets." },
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyTypeObject glog_WizardType = {
    PyVarObject_HEAD_INIT(NULL, 0)
        "glog.Wizard",             /* tp_name */
    sizeof(glog_Wizard),             /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor) wizard_dealloc, /* tp_dealloc */
    0,                         /* tp_print */
    0,                         /* tp_getattr */
    0,                         /* tp_setattr */
    0,                         /* tp_reserved */
    0,                         /* tp_repr */
    0,                         /* tp_as_number */
    0,                         /* tp_as_sequence */
    0,                         /* tp_as_mapping */
    0,                         /* tp_hash  */
    0,                         /* tp_call */
    0,                         /* tp_str */
    0,                         /* tp_getattro */
    0,                         /* tp_setattro */
    0,                         /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT |
        Py_TPFLAGS_BASETYPE,   /* tp_flags */
    "Wizard",           /* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    Wizard_methods,             /* tp_methods */
    0,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)wizard_init,      /* tp_init */
    0,                         /* tp_alloc */
    wizard_new,                 /* tp_new */
};

static PyObject * wizard_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
    glog_Wizard *self;
    self = (glog_Wizard*)type->tp_alloc(type, 0);
    return (PyObject *)self;
}

static int wizard_init(glog_Wizard *self, PyObject *args, PyObject *kwds) {
    return 0;
}

static void wizard_dealloc(glog_Wizard* self) {
    Py_TYPE(self)->tp_free((PyObject*)self);
}

extern PyTypeObject glog_ProgramType;
static PyObject *wizard_rewrite_program(PyObject *self, PyObject *args) {
    PyObject *program = NULL;
    const char *query = NULL;
    if (PyArg_ParseTuple(args, "Os", &program, &query)) {
        if (strcmp(program->ob_type->tp_name, "glog.Program") == 0) {
            glog_Program *p = (glog_Program*)program;
            auto program = p->program;
            //Parse query
            auto sQuery = std::string(query);
            Dictionary dv;
            auto lQuery = program->parseLiteral(sQuery, dv);

            Wizard w;
            auto adornedProgram = w.getAdornedProgram(lQuery, *program.get());

            std::pair<PredId_t, PredId_t> ioPredIDs;
            auto newProgram = w.doMagic(lQuery, adornedProgram, ioPredIDs);

            //Return a new program
            auto arglist = Py_BuildValue("(O)", (PyObject*)p->e);
            PyObject *obj = PyObject_CallObject((PyObject *) &glog_ProgramType,
                    arglist);
            Py_DECREF(arglist);
            glog_Program *newP = (glog_Program*)obj;
            newP->program = newProgram;
            newP->ioPredIDs = ioPredIDs;
            return obj;
        }
    }
    Py_INCREF(Py_None);
    return Py_None;

}
