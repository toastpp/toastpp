/*
 *  Copyright 2008-2009 NVIDIA Corporation
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */



#pragma once
#include <string.h>

char * get_arg(int argc, char ** argv, const char * key)
{
    char * val = NULL;

    size_t keylen = strlen(key);

    for(int i = 1; i < argc; i++){
        char * token = argv[i];
        
        if(strncmp(token, "--", 2) != 0)
            continue;
        token += 2;

        if(strncmp(token, key, keylen) != 0)
            continue;
        token += keylen;

        val = argv[i];
    }
    
    return val;
}

char * get_argval(int argc, char ** argv, const char * key)
{
    char * val = NULL;

    size_t keylen = strlen(key);

    for(int i = 1; i < argc; i++){
        char * token = argv[i];
        
        if(strncmp(token, "--", 2) != 0)
            continue;
        token += 2;

        if(strncmp(token, key, keylen) != 0)
            continue;
        token += keylen;

        if(strncmp(token, "=", 1) != 0)
            continue;
        token += 1;

        val = token;
    }
    
    return val;
}

