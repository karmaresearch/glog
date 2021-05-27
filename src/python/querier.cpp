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

/*** Methods ***/
static PyObject * querier_new(PyTypeObject *type, PyObject *args, PyObject *kwds);
static int querier_init(glog_Querier *self, PyObject *args, PyObject *kwds);
static void querier_dealloc(glog_Querier* self);
static PyObject* querier_get_derivation_tree(PyObject* self, PyObject *args);

static PyMethodDef Querier_methods[] = {
    {"get_derivation_tree", querier_get_derivation_tree, METH_VARARGS, "Get derivation tree of a fact." },

    {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyTypeObject glog_QuerierType = {
    PyVarObject_HEAD_INIT(NULL, 0)
        "glog.Querier",             /* tp_name */
    sizeof(glog_Querier),             /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor) querier_dealloc, /* tp_dealloc */
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
    "Querier",           /* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    Querier_methods,             /* tp_methods */
    0,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)querier_init,      /* tp_init */
    0,                         /* tp_alloc */
    querier_new,                 /* tp_new */
};

static PyObject * querier_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
    glog_Querier *self;
    self = (glog_Querier*)type->tp_alloc(type, 0);
    self->g = NULL;
    return (PyObject *)self;
}

static int querier_init(glog_Querier *self, PyObject *args, PyObject *kwds) {
    PyObject *arg = NULL;
    if (!PyArg_ParseTuple(args, "O", &arg))
        return -1;
    if (arg != NULL) {
        if (strcmp(arg->ob_type->tp_name, "glog.TG") != 0)
            return -1;
        Py_INCREF(arg);
        self->g = (glog_TG*)arg;
        self->q = self->g->g->getQuerier();
    }
    return 0;
}

static void querier_dealloc(glog_Querier* self) {
    if (self->g != NULL) {
        Py_DECREF(self->g);
    }
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject* querier_get_derivation_tree(PyObject* self, PyObject *args) {
    size_t nodeId;
    size_t factId;
    if (!PyArg_ParseTuple(args, "ll", &nodeId, &factId)) {
        Py_INCREF(Py_None);
        return Py_None;
    }
    auto out = ((glog_Querier*)self)->q->getDerivationTree(nodeId, factId);
    std::stringstream ssOut;
    JSON::write(ssOut, out);
    std::string sOut = ssOut.str();
    //Return JSON object
    return PyUnicode_FromString(sOut.c_str());
}
