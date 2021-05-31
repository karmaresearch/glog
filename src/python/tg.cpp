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
static PyObject * tg_new(PyTypeObject *type, PyObject *args, PyObject *kwds);
static int tg_init(glog_TG *self, PyObject *args, PyObject *kwds);
static void tg_dealloc(glog_TG* self);
static PyObject* tg_add_node(PyObject* self, PyObject *args);
static PyObject* tg_get_n_nodes(PyObject* self, PyObject *args);
static PyObject* tg_get_n_edges(PyObject* self, PyObject *args);
static PyObject* tg_get_n_facts(PyObject* self, PyObject *args);
static PyObject* tg_get_node_size(PyObject* self, PyObject *args);

static PyMethodDef TG_methods[] = {
    {"get_n_nodes", tg_get_n_nodes, METH_VARARGS, "Get n. nodes in the TG." },
    {"get_n_edges", tg_get_n_edges, METH_VARARGS, "Get n. edges in the TG." },
    {"get_n_facts", tg_get_n_facts, METH_VARARGS, "Get n. facts in the TG." },
    {"get_node_size", tg_get_node_size, METH_VARARGS, "Get number of facts stored in a node." },
    {"add_node", tg_add_node, METH_VARARGS, "Add a node with some provided facts." },
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyTypeObject glog_TGType = {
    PyVarObject_HEAD_INIT(NULL, 0)
        "glog.TG",             /* tp_name */
    sizeof(glog_TG),             /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor) tg_dealloc, /* tp_dealloc */
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
    "TG",           /* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    TG_methods,             /* tp_methods */
    0,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)tg_init,      /* tp_init */
    0,                         /* tp_alloc */
    tg_new,                 /* tp_new */
};

static PyObject * tg_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
    glog_TG *self;
    self = (glog_TG*)type->tp_alloc(type, 0);
    self->reasoner = NULL;
    self->g = NULL;
    return (PyObject *)self;
}

static int tg_init(glog_TG *self, PyObject *args, PyObject *kwds) {
    PyObject *arg = NULL;
    if (!PyArg_ParseTuple(args, "O", &arg))
        return -1;
    if (arg != NULL) {
        if (strcmp(arg->ob_type->tp_name, "glog.Reasoner") != 0)
            return -1;
        Py_INCREF(arg);
        self->reasoner = (glog_Reasoner*)arg;
        auto &g = self->reasoner->sn->getGBGraph();
        self->g = &g;
    }
    return 0;
}

static void tg_dealloc(glog_TG* self) {
    if (self->reasoner != NULL) {
        Py_DECREF(self->reasoner);
    }
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject* tg_add_node(PyObject* self, PyObject *args) {
    PredId_t predid;
    size_t step;
    PyObject *facts;
    std::vector<std::vector<std::string>> parsedFacts;
    if (PyArg_ParseTuple(args, "ILO", &predid, &step, &facts)) {
        PyObject *iter = PyObject_GetIter(facts);
        while (true) {
            PyObject *next = PyIter_Next(iter);
            if (!next) {
                break;
            }
            if (PyTuple_Check(next)) {
                parsedFacts.emplace_back();
                auto card = PyTuple_Size(next);
                for(size_t i = 0; i < card; ++i) {
                    auto el = PyTuple_GetItem(next, i);
                    const char *s = PyBytes_AsString(el);
                    parsedFacts.back().push_back(std::string(s));
                }
            }
            Py_DECREF(next);
        }
    }
    if (!parsedFacts.empty()) {
        glog_TG *s = (glog_TG*)self;
        s->g->addNode(predid, step, parsedFacts);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject* tg_get_n_nodes(PyObject* self, PyObject *args) {
    glog_TG *s = (glog_TG*)self;
    return PyLong_FromLong(s->g->getNNodes());
}

static PyObject* tg_get_n_facts(PyObject* self, PyObject *args) {
    glog_TG *s = (glog_TG*)self;
    return PyLong_FromLong(s->g->getNNodes());
}

static PyObject* tg_get_n_edges(PyObject* self, PyObject *args) {
    glog_TG *s = (glog_TG*)self;
    return PyLong_FromLong(s->g->getNEdges());
}

static PyObject* tg_get_node_size(PyObject* self, PyObject *args) {
    glog_TG *s = (glog_TG*)self;
    size_t nodeId;
    if (!PyArg_ParseTuple(args, "l", &nodeId)) {
        Py_INCREF(Py_None);
        return Py_None;
    }
    return PyLong_FromLong(s->g->getNodeSize(nodeId));
}
