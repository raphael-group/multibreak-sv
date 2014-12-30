import numpy as np
cimport numpy as np

ctypedef np.int8_t dtype_ti
ctypedef np.float_t dtype_tf


def outputDelInsCython (bytes alignedQuery, bytes alignedTarget, np.ndarray[dtype_tf, ndim=2] t_log2, np.ndarray[dtype_tf, ndim=2] e_log2):
    cdef int n = len(alignedQuery)
    cdef np.ndarray[dtype_tf, ndim=2] v = np.zeros((3,n+2), dtype=float)
    cdef np.ndarray[dtype_ti, ndim=2] p = np.zeros((3,n+2), dtype='int8')

    cdef int i = 0
    cdef int j = 0
    cdef int x = 1
    maxval = -1
    argmax = 0
    cdef int prevState = -1
    cdef int currState = 2
    cdef np.ndarray[dtype_tf, ndim=2] temp =  np.zeros((3,3), dtype=float)
    cdef np.ndarray[dtype_ti, ndim=1]  states = np.zeros(n, dtype='int8')
    cdef int seen2Tran = 0

    # set all pointers to point to normal case in second row
    v[2][0] = 0.0
    v[0][0] = -1000.0
    v[1][0] = -1000.0

    for i in range(n):
        x = 1
        if alignedQuery[i] == '-':
            x = 2
        if alignedTarget[i] == '-':
            x = 0
	
        for j in range (2):
            maxval = v[0,i] + t_log2[j,0]
            maxarg = 0
            if maxval <= v[1,i] + t_log2[j,1]:
                maxval = v[1,i] + t_log2[j,1]
                maxarg = 1
            if maxval <= v[2,i] + t_log2[j,2]:
                maxval = v[2,i] + t_log2[j,2]
                maxarg = 2
            v[j,i+1] = e_log2[j,x]+maxval
            p[j,i+1] = maxarg
        # below is an unnecessaery speed up to eliminate checking
        # the normal transition case (basically reduces an if call 
        # in the loop
        j = 2
        maxval = v[0,i] + t_log2[j,0]
        maxarg = 0
        if maxval <= v[1,i] + t_log2[j,1]:
            maxval = v[1,i] + t_log2[j,1]
            maxarg = 1
        if maxval <= v[2,i] + t_log2[j,2]:
            maxval = v[2,i] + t_log2[j,2]
            maxarg = 2
        v[j,i+1] = e_log2[j,x]+maxval
        p[j,i+1] = maxarg
        if maxarg != 2:
            seen2Tran = 1
    
    if seen2Tran == 0:
        return [0]

    states[n-1] = 2
    for i in range (n, 1, -1):
        prevState = p[currState][i]
        states[i-2] = prevState
        currState = prevState

    return states

def printSummaryFast (bytes qid, bytes tid, int t_s, int t_e, bytes strand, np.ndarray[dtype_ti, ndim=1] states, bytes alignedQuery, bytes alignedTarget, out):
    cdef int n = len(states)
    cdef int normal = 1
    cdef int event_start = -1
    cdef int event_pos = -1
    cdef char event = 'N'
    cdef int pos = t_s
    cdef int event_size = 0
    cdef int i = 0
    cdef int qpos = 0
    cdef int qstart = 0
    cdef int qinc = 1
    if strand[0] == '-':
        qpos = n-1
        qinc = -1
    cdef int qend = qstart
    for i in range (n):
        if states[i] != 2:
            if normal:
                normal = 0
                event_start = pos
                event_pos = pos
                qstart = qpos
                qend = qpos
                eventsize = 1
                if states[i]:
                    event = 'D'
                else:
                    event = 'I'
            else:
                event_pos = pos
                event_size += 1
                qend = qpos
        else:
            if not(normal):
                qdist = abs(qend-qstart)
                out.write("%s %s %i %i %i %i %s %c %i %i %i %i\n" %( qid, tid, t_s, t_e, event_start, event_pos+1, strand, event, event_size, qdist, qstart, qend ))
            event_size = 0
            normal = 1
        if alignedTarget[i] != "-":
            pos += 1
        if alignedQuery[i] != "-":
            qpos += qinc
        
