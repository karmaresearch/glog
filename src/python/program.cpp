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

#include <python/glog.h>
#include <kognac/utils.h>

/*** Methods ***/
static PyObject * program_new(PyTypeObject *type, PyObject *args, PyObject *kwds);
static int program_init(glog_Program *self, PyObject *args, PyObject *kwds);
static void program_dealloc(glog_Program* self);
static PyObject* program_load_from_file(PyObject* self, PyObject *args);
static PyObject* program_get_n_rules(PyObject* self, PyObject *args);
static PyObject* program_get_rule(PyObject* self, PyObject *args);
static PyObject* program_get_predicate_name(PyObject* self, PyObject *args);

static PyMethodDef Program_methods[] = {
    {"load_from_file", program_load_from_file, METH_VARARGS, "Load rules from file." },
    {"get_n_rules", program_get_n_rules, METH_VARARGS, "Return n rules." },
    {"get_rule", program_get_rule, METH_VARARGS, "Get rules." },
    {"get_rule", program_get_rule, METH_VARARGS, "Get rules." },
    {"get_predicate_name", program_get_predicate_name, METH_VARARGS, "Get name predicate." },
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyTypeObject glog_ProgramType = {
    PyVarObject_HEAD_INIT(NULL, 0)
        "glog.Program",             /* tp_name */
    sizeof(glog_Program),             /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor) program_dealloc, /* tp_dealloc */
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
    "Program",           /* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    Program_methods,             /* tp_methods */
    0,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)program_init,      /* tp_init */
    0,                         /* tp_alloc */
    program_new,                 /* tp_new */
};

static PyObject * program_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
    glog_Program *self;
    self = (glog_Program*)type->tp_alloc(type, 0);
    self->e = NULL;
    self->program = NULL;
    return (PyObject *)self;
}

static int program_init(glog_Program *self, PyObject *args, PyObject *kwds) {
    PyObject *arg = NULL;
    if (!PyArg_ParseTuple(args, "O", &arg))
        return -1;
    if (arg != NULL) {
        if (strcmp(arg->ob_type->tp_name, "glog.EDBLayer") != 0)
            return -1;
        Py_INCREF(arg);
        self->e = (glog_EDBLayer*)arg;
        self->program = std::shared_ptr<Program>(new Program(self->e->e));
    }
    return 0;
}

static void program_dealloc(glog_Program* self) {
    if (self->program) {
        Py_DECREF(self->e);
    }
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject* program_load_from_file(PyObject *self, PyObject *args) {
    const char *path = NULL;
    if (PyArg_ParseTuple(args, "|s", &path)) {
        ((glog_Program*)self)->program->readFromFile(std::string(path));
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject* program_get_n_rules(PyObject* self, PyObject *args) {
    auto nrules = ((glog_Program*)self)->program->getNRules();
    return PyLong_FromLong(nrules);
}

static PyObject* program_get_rule(PyObject* self, PyObject *args) {
    size_t ruleIdx = 0;
    if (PyArg_ParseTuple(args, "i", &ruleIdx)) {
        //Get the rule
        auto rule =  ((glog_Program*)self)->program->getRule(ruleIdx);
        std::string sRule = rule.tostring(((glog_Program*)self)->program.get(),
                ((glog_Program*)self)->e->e);
        return PyUnicode_FromStringAndSize(sRule.c_str(), sRule.size());
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject* program_get_predicate_name(PyObject* self, PyObject *args) {
    size_t predId = 0;
    if (!PyArg_ParseTuple(args, "l", &predId)) {
        Py_INCREF(Py_None);
        return Py_None;
    }
    auto predName = ((glog_Program*)self)->program->getPredicateName(predId);
    return PyUnicode_FromString(predName.c_str());
}
