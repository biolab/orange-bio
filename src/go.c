#include<Python.h>
#include<structmember.h>
#include<string.h>
#include<stdio.h>

#define RAISE_TYPE_ERROR(error) PyErr_SetString(PyExc_TypeError, error); Py_INCREF(PyExc_TypeError);
#ifndef MAX
#define MAX(a,b) (((a)<(b))?(b):(a))
#endif

struct _TermNode;
typedef struct _TermNode{
	char* id;
	struct _TermNode* next;
} TermNode;

TermNode* makeTermNode(char* id){
	TermNode* tmp=malloc(sizeof(TermNode));
	tmp->id=id;
	tmp->next=NULL;
	return tmp;
}

typedef struct{
	TermNode* head;
	TermNode* tail;
} TermList;

TermList* makeTermList(){
	TermList *list=malloc(sizeof(TermList));
	memset(list, 0, sizeof(TermList));
	return list;
}

void addTermNode(TermList* list, TermNode* node){
	TermNode* tmp=list->head;
	list->head=node;
	node->next=tmp;
	if(!list->tail)
		list->tail=node;
}

void clearTermList(TermList* list){
	TermNode* node=list->head;
	while(node){
		TermNode* tmp=node;
		node=node->next;
		free(tmp);
	}
	list->head=NULL;
	list->tail=NULL;
}

struct _GeneEvidenceNode;
typedef struct _GeneEvidenceNode{
	char *geneName;
	int evidence;
	struct _GeneEvidenceNode* next;
} GeneEvidenceNode;

GeneEvidenceNode* makeGeneEvidenceNode(char* name, int evidence){
	GeneEvidenceNode* tmp=malloc(sizeof(GeneEvidenceNode));
	tmp->geneName=name;
	tmp->evidence=evidence;
	tmp->next=NULL;
	return tmp;
}

typedef struct {
	GeneEvidenceNode* head;
	GeneEvidenceNode* tail;
} GeneEvidenceList;

GeneEvidenceList* makeGeneEvidenceList(){
	GeneEvidenceList *list=malloc(sizeof(GeneEvidenceList));
	memset(list, 0, sizeof(GeneEvidenceList));
	return list;
}

void addGeneEvidenceNode(GeneEvidenceList* list, GeneEvidenceNode* node){
	GeneEvidenceNode* tmp=list->head;
	list->head=node;
	node->next=tmp;
	if(!list->tail)
		list->tail=node;
}

void clearGeneEvidenceList(GeneEvidenceList* list){
	GeneEvidenceNode* node=list->head;
	while(node){
		GeneEvidenceNode* tmp=node;
		node=node->next;
		free(tmp);
	}

	list->head=NULL;
	list->tail=NULL;
}

struct _hash_entry;
typedef struct _hash_entry{
	char *key;
	void *ptr;
	struct _hash_entry* next;
}hash_entry;

typedef struct{
	int buckets;
	hash_entry * entries;
}hash_table;

hash_entry* makeHashEntry(char* key, void* ptr){
	hash_entry* tmp=malloc(sizeof(hash_entry));
	tmp->key=strdup(key);
	tmp->ptr=ptr;
	tmp->next=NULL;
	return tmp;
}

unsigned int hash_fun(char* key, int maxVal){
	unsigned int sum=0;
	maxVal=MAX(maxVal, 1);
	while(*key){
		sum=sum*131+ *key;
		key++;
	}
	return sum%maxVal;
}

void* hash_add(hash_table* hash, char* key, void* ptr){
	unsigned int hashVal=hash_fun(key, hash->buckets);
	hash_entry* entry=&hash->entries[hashVal];
	if(entry->key){	//conflict
		while(entry->next){
			/*if the key allready in the hash table replace the ptr and 
			return the old ptr*/
			if(!strcmp(entry->key, key)){ 
				void *tmp=entry->ptr;
				entry->ptr=ptr;
				printf("replacing hash entry %s\n", key);
				return tmp;
			}
			entry=entry->next;
		}
		entry->next=makeHashEntry(key, ptr);
	}
	else {
		entry->key=strdup(key);
		entry->ptr=ptr;
		entry->next=NULL;
	}
	return NULL;
}

#define HASH_MISS (void*)-1

void *hash_get(hash_table* hash, char* key){
	int hashVal=hash_fun(key, hash->buckets);
	hash_entry* entry=&hash->entries[hashVal];
	//printf("hash of %s:%i\n", key, hashVal);
	while(entry && entry->key)
		if(strcmp(key, entry->key)==0)
			return entry->ptr;
		else{
			//printf("%s", entry->key);
			entry=entry->next;
		}
	//printf("HASH MISS: %s\n",key);
	return HASH_MISS;
}

char hash_has_key(hash_table* hash, char* key){
	return hash_get(hash, key)!=HASH_MISS;
}

hash_table* makeHashTable(int buckets){
	hash_table* hash=malloc(sizeof(hash_table));
	hash->buckets=MAX(buckets,1);
	hash->entries=malloc(sizeof(hash_entry)*buckets);
	memset(hash->entries, 0, sizeof(hash_entry)*buckets);
	return hash;
}

void freeHashTable(hash_table* table){
	hash_entry* entry=NULL;
	hash_entry* tmp=NULL;
	int i;
	for(i=0;i<table->buckets;i++){
		entry=&table->entries[i];
		if(entry->key)
			free(entry->key);
		entry=entry->next;
		while(entry){
			free(entry->key);
			tmp=entry;
			entry=entry->next;
			free(tmp);
		}
	}
	free(table->entries);
	free(table);
}

typedef struct{
	char* goID;
	GeneEvidenceList mappedGenes;
	TermList parents;
	TermList children;
	int numRef;
	int code;
	char visited;
} GOTerm;

typedef struct{
	char* name;
	char* goID;
	int aspect;
	int evidence;
} AnnRecord;

typedef struct{
	PyObject_HEAD
	PyObject* aliasMapper;
	hash_table* hash;
	AnnRecord* annotation;
	int numGenes;
} Annotation;

static PyMemberDef Annotation_members[]={
	{"aliasMapper", T_OBJECT_EX, offsetof(Annotation,aliasMapper), 0, ""},
	{NULL}
};
typedef struct{
	PyObject_HEAD
	PyObject* aliasMapper;
	int aspect;
	hash_table* hash;
	GOTerm *terms;
} Ontology;

static PyMemberDef Ontology_members[]={
	{"aliasMapper", T_OBJECT_EX, offsetof(Ontology, aliasMapper), 0, ""},
	{NULL}
};

void prepareGOTerms(Ontology* ontology);

static void Annotation_dealloc(Annotation* self){
	AnnRecord* rec=self->annotation;
	//printf("Annotation.dealloc\n");
	while(rec->name){
		free(rec->name);
		free(rec->goID);
		rec++;
	}
	free(self->annotation);
	freeHashTable(self->hash);
	Py_XDECREF(self->aliasMapper);
	self->ob_type->tp_free(self);
}

static void Ontology_dealloc(Ontology* self){
	GOTerm* term=NULL;
	TermNode* node=NULL;
	//printf("Ontology.dealloc\n");
	prepareGOTerms(self);
	term=self->terms;
	while(term->goID){
		/* free the parent list of the term*/
		node=term->parents.head;
		while(node){
			free(node->id);
			node=node->next;
		}
		clearTermList(&term->parents);
		clearTermList(&term->children);

		free(term->goID);
		term++;
	}
	freeHashTable(self->hash);
	free(self->terms);
	Py_XDECREF(self->aliasMapper);
	self->ob_type->tp_free(self);
}

static PyTypeObject go_AnnotationType={
	PyObject_HEAD_INIT(NULL)
	0,                         /*ob_size*/
    "go.Annotation",           /*tp_name*/
    sizeof(Annotation),        /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)Annotation_dealloc,  /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT,        /*tp_flags*/
    "Holds the annotation records",      /*tp_doc*/
	0,		               /* tp_traverse */
    0,		               /* tp_clear */
    0,		               /* tp_richcompare */
    0,		               /* tp_weaklistoffset */
    0,		               /* tp_iter */
    0,		               /* tp_iternext */
    0,                     /* tp_methods */
    Annotation_members,             /* tp_members */

};

static PyTypeObject go_OntologyType={
	PyObject_HEAD_INIT(NULL)
	0,                         /*ob_size*/
    "go.Ontology",               /*tp_name*/
    sizeof(Ontology),            /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)Ontology_dealloc,      /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT,        /*tp_flags*/
    "Holds the go terms",      /*tp_doc*/
	0,		               /* tp_traverse */
    0,		               /* tp_clear */
    0,		               /* tp_richcompare */
    0,		               /* tp_weaklistoffset */
    0,		               /* tp_iter */
    0,		               /* tp_iternext */
    0,                     /* tp_methods */
    Ontology_members,             /* tp_members */
};

GOTerm* getGOTerm(Ontology* ontology, char* id){
	GOTerm* term=hash_get(ontology->hash, id);
	if(term==HASH_MISS && PyMapping_Check(ontology->aliasMapper) && PyMapping_HasKeyString(ontology->aliasMapper, id)){
		PyObject* string=PyMapping_GetItemString(ontology->aliasMapper, id);
		id=PyString_AsString(string);
		term=hash_get(ontology->hash, id);
		Py_DECREF(string);
	}
	/*if(term==HASH_MISS)
		printf("HASH MISS: %s\n",id);*/

	return term;
}

AnnRecord* getAnnRecord(Annotation* annotation, char* gene){
	AnnRecord* rec=hash_get(annotation->hash, gene);
	if(rec==HASH_MISS && PyMapping_Check(annotation->aliasMapper) && PyMapping_HasKeyString(annotation->aliasMapper, gene)){
		PyObject* string=PyMapping_GetItemString(annotation->aliasMapper, gene);
		gene=PyString_AsString(string);
		rec=hash_get(annotation->hash, gene);
		Py_DECREF(string);
	}
	/*if(rec==HASH_MISS)
		printf("HASH MISS: %s\n",gene);*/
	return rec;
}

GOTerm makeGOTerm(PyObject* tuple){
	GOTerm term;
	memset(&term, 0, sizeof(GOTerm));
	if(tuple && PyTuple_Check(tuple)){
		char* id=strdup(PyString_AsString(PyTuple_GetItem(tuple, 0)));
		PyObject* parents=PyTuple_GetItem(tuple, 1);
		int nParents=0;
		TermList parentList;
		parentList.head=NULL;
		parentList.tail=NULL;
		if(PyList_Check(parents) && (nParents=PyList_Size(parents))){
			PyObject* iter=PyObject_GetIter(parents);
			PyObject* item=NULL;
			while((item=PyIter_Next(iter))){
				addTermNode(&parentList, makeTermNode(strdup(PyString_AsString(item))));
				Py_DECREF(item);
			}
			Py_DECREF(iter);
		}
		term.goID=id;
		term.parents=parentList;
		term.children.head=NULL;
		term.children.tail=NULL;
		return term;
	}
	else{
		printf("Could not create a proper GOTerm");
		term.goID="Uninitialized GO Term";
		term.parents.head=NULL;
		term.parents.tail=NULL;
		term.children.head=NULL;
		term.children.tail=NULL;
		return term;
	}
}	

AnnRecord makeAnnRecord(PyObject* tuple){
	AnnRecord ann;
	memset(&ann, 0, sizeof(AnnRecord));
	if(tuple && PyTuple_Check(tuple)){
		ann.name=strdup(PyString_AsString(PyTuple_GetItem(tuple, 0)));
		ann.goID=strdup(PyString_AsString(PyTuple_GetItem(tuple, 1)));
		ann.evidence=PyInt_AsLong(PyTuple_GetItem(tuple, 2));
		ann.aspect=PyInt_AsLong(PyTuple_GetItem(tuple,3));
	} else{
		printf("Not a tuple. Failed to build a proper AnnRecord\n");
		ann.name="Uninitialized annotation record";
		ann.goID=NULL;
		ann.evidence=0;
		ann.aspect=0;
	}
	return ann;
}

PyObject* parseGOTerms(PyObject* self, PyObject* args){
	PyObject* terms=NULL;
	PyObject* iter=NULL;
	PyObject* item=NULL;
	PyObject* aliasMapper=Py_None;
	GOTerm* goTerms=NULL;
	GOTerm* ptr=NULL;
	hash_table* hash=NULL;
	Ontology *base=NULL;
	int nTerms=0;
	int maxTerms=0;
	if(!PyArg_ParseTuple(args, "O|O:parseGOTerms", &terms, &aliasMapper)){
		return NULL;
	}
	if(!PyList_Check(terms)){
		RAISE_TYPE_ERROR("Wrong argument type. Need a list");
		return NULL;
	}
	nTerms=PyList_Size(terms);
	maxTerms=nTerms+1;
	goTerms=malloc(sizeof(GOTerm)*maxTerms);
	memset(goTerms, 0, sizeof(GOTerm)*maxTerms);
	hash=makeHashTable(nTerms*3/2);
	iter=PyObject_GetIter(terms);
	ptr=goTerms;
	while((item=PyIter_Next(iter))){
		(*ptr)=makeGOTerm(item);
		hash_add(hash, ptr->goID, ptr);
		Py_DECREF(item);
		ptr++;
	}
	Py_DECREF(iter);
	base=(Ontology*)go_OntologyType.tp_new(&go_OntologyType, NULL, NULL);
	base->terms=goTerms;
	base->hash=hash;
	base->aliasMapper=aliasMapper;
	Py_INCREF(aliasMapper);
	ptr=goTerms;
	while(ptr->goID){
		TermNode* node=ptr->parents.head;
		while(node){
			GOTerm* term=getGOTerm(base, node->id);
			if(term!=HASH_MISS)
				addTermNode(&term->children, makeTermNode(ptr->goID));
			node=node->next;
		}
		ptr++;
	}
	return (PyObject*)base;
}
	
PyObject* parseAnnotation(PyObject* self, PyObject* args){
	PyObject* ann=NULL;
	PyObject* iter=NULL;
	PyObject* item=NULL;
	PyObject* aliasMapper=Py_None;
	AnnRecord* annotation=NULL;
	AnnRecord* ptr=NULL;
	Annotation* annot=NULL;
	hash_table *hash=NULL;

	char* lastName="I'm Bender baby! Please insert liquor!";
	int nAnn=0;
	int annSize=0;
	int numGenes=0;
	if(!PyArg_ParseTuple(args, "O|O:parseAnnotation", &ann, &aliasMapper)){
		return NULL;
	}
	if(!PyList_Check(ann)){
		RAISE_TYPE_ERROR("Need a list");
		return NULL;
	}
	iter=PyObject_GetIter(ann);
	nAnn=PyList_Size(ann);
	annSize=nAnn+1;
	annotation=malloc(sizeof(AnnRecord)*annSize);
	memset(annotation, 0, sizeof(AnnRecord)*annSize);
	hash=makeHashTable(annSize);
	ptr=annotation;
	while((item=PyIter_Next(iter))){
		*ptr=makeAnnRecord(item);
		if(strcmp(ptr->name, lastName)){
			hash_add(hash, ptr->name, ptr);
			lastName=ptr->name;
			numGenes++;
		}
		ptr++;
		Py_DECREF(item);
	}
	Py_DECREF(iter);
	annot=(Annotation*)go_AnnotationType.tp_new(&go_AnnotationType, NULL, NULL);
	annot->annotation=annotation;
	annot->hash=hash;
	annot->aliasMapper=aliasMapper;
	annot->numGenes=numGenes;
	Py_INCREF(aliasMapper);
	return (PyObject*)annot; 
}

void prepareGOTerms(Ontology* base){
	GOTerm* term=base->terms;
	while(term->goID){
		term->visited=0;
		term->code=0;
		term->numRef=0;
		clearGeneEvidenceList(&term->mappedGenes);
		term++;
	}
}

void markTerms(GOTerm* term, Ontology* ontology, char* geneName, int evidence, char recursive, TermList* selected){
	/*check if the gene allready on the selected list*/
    TermNode *tnode=selected->head;
	GeneEvidenceNode* gnode=term->mappedGenes.head;
    while(tnode && strcmp(tnode->id, term->goID))
		tnode=tnode->next;
    /*if not add it*/
	if(!tnode){
		//printf("selected node %s\n", term->goID);
		addTermNode(selected, makeTermNode(term->goID));
	}
		
	/*add the gene/evidence to the term */
	
	while(gnode && strcmp(gnode->geneName, geneName))
		gnode=gnode->next;
	if(gnode)
		gnode->evidence|=evidence;
	else
		/*add the gene name to term*/
		addGeneEvidenceNode(&term->mappedGenes, makeGeneEvidenceNode(geneName, evidence));
	//printf("marking term %s\n", term->goID);

	if(recursive){
		tnode=term->parents.head;
		while(tnode){
			if((term=getGOTerm(ontology, tnode->id))!=HASH_MISS)
				markTerms(term, ontology, geneName, evidence, recursive, selected);
			tnode=tnode->next;
		}
	}
}

#define DIRECT_MAPPING 1
#define INDIRECT_MAPPING 2

void clearVisited(Ontology* go){
	GOTerm* term=go->terms;
	while(term->goID){
		term->visited=0;
		term++;
	}
}

/*recursivly mark ALL parents of a term to map indirectly, including
 the ones that are allready marked as direct  */
void _slimMapping2(GOTerm* term, Ontology* ontology, Ontology* slimOntology, TermList* list){
	GOTerm* parent=NULL;
	TermNode* parentNode=NULL;
	parentNode=term->parents.head;
	while(parentNode){
		if((parent=getGOTerm(ontology, parentNode->id))==HASH_MISS){
			parentNode=parentNode->next;
			continue;
		}
		parent->visited=INDIRECT_MAPPING;
		_slimMapping2(parent, ontology, slimOntology, list);
		parentNode=parentNode->next;
	}
}

void _slimMapping1(GOTerm* term, Ontology* ontology, Ontology* slimOntology, TermList* list){
	GOTerm* parent=NULL;
	TermNode* parentNode=NULL;
	parentNode=term->parents.head;
	while(parentNode){
		if((parent=getGOTerm(ontology, parentNode->id))==HASH_MISS){
			parentNode=parentNode->next;
			continue;
		}
		
		if(getGOTerm(slimOntology, parent->goID)!=HASH_MISS && parent->visited==0){
			/*if the parent in the slims subset and not yet visited
			mark it and add it to the list and recursivly mark ALL parents of 
			a term to map indirectly, including the ones that are allready marked as direct*/
			parent->visited=DIRECT_MAPPING;
			addTermNode(list, makeTermNode(parent->goID));
			_slimMapping2(parent, ontology, slimOntology, list);
		}else if(parent->visited==0)
			/*else if not yet visited, NOTE:do nothing if term allready visited*/
			_slimMapping1(parent, ontology, slimOntology, list);
		parentNode=parentNode->next;
	}
}


void slimMapping(GOTerm* term, Ontology* ontology, Ontology* slimOntology, TermList* list){
	clearVisited(ontology);
	if(getGOTerm(slimOntology, term->goID)!=HASH_MISS){
		/*the term is allready in the slim subset*/
		term->visited=DIRECT_MAPPING;
		addTermNode(list, makeTermNode(term->goID));
		return;
	}
	//printf("mapping slim: %s\n", goId);
	_slimMapping1(term, ontology, slimOntology, list);
}

TermList* findTerms(char** geneName, int evidence, int aspect, char slimsOnly, char recursive, Annotation* annotation, Ontology* ontology, Ontology* slimOntology, PyObject* callback, int step, int start){
	TermList* list=malloc(sizeof(TermList));
	GOTerm* term=NULL;
	int count=0;
	list->head=NULL;
	list->tail=NULL;

	while(*geneName){
		AnnRecord* ann=NULL;
		char *mappedGeneName=NULL;
		count++;
		if(callback && count%step==0 && start<=100)
			PyObject_CallFunction(callback, "i", start++);

		if((ann=getAnnRecord(annotation, *geneName))==HASH_MISS){
			printf("unknown geneName %s\n", *geneName);
			geneName++;
			continue;
		}
		
		mappedGeneName=ann->name;
		while(ann->name && !strcmp(ann->name, mappedGeneName)){
			if((ann->aspect & aspect) && (evidence & ann->evidence)){
				TermList termList;
				TermNode* node=NULL;
				if((term=getGOTerm(ontology, ann->goID))==HASH_MISS){
					ann++;
					continue;
				}
				memset(&termList, 0, sizeof(TermList));
				if(slimsOnly && slimOntology){
					slimMapping(term, ontology, slimOntology, &termList);
					node=termList.head;
					while(node){
						if((term=getGOTerm(ontology, node->id))!=HASH_MISS && term->visited==DIRECT_MAPPING)
							markTerms(term, ontology, *geneName, ann->evidence, recursive, list);
						node=node->next;
					}
				}else
					markTerms(term, ontology, *geneName, ann->evidence, recursive, list);
				clearTermList(&termList);
			}
			ann++;
		}
		geneName++;
	}
	return list;
}

void addReference(GOTerm* term, Ontology* ontology, int code){
	TermNode* node=NULL;
	if(term->code==code)
		return;
	term->numRef++;
	term->code=code;
	//printf("adding ref %s\n", term->goID);

	node=term->parents.head;
	while(node){
		if((term=getGOTerm(ontology, node->id))!=HASH_MISS)
			addReference(term, ontology, code);
		node=node->next;
	}
}

void mapReferenceGenes(char** geneName, int evidence, int aspect, Annotation* annotation, Ontology* ontology, PyObject* callback, int step, int start){
	GOTerm* term=NULL;
	int markCode=1;
	int count=0;
	while(*geneName){
		AnnRecord* ann=NULL;
		char *mappedGeneName=NULL;
		count++;
		if(callback && count%step==0)
			PyObject_CallFunction(callback, "i",start++);

		if((ann=getAnnRecord(annotation, *geneName))==HASH_MISS){
			printf("unknown geneName %s\n", *geneName);
			geneName++;
			continue;
		}
		mappedGeneName=ann->name;
		//printf("%s, %s\n", *geneName, ann->name);
		while(ann->name && strcmp(ann->name, mappedGeneName)==0){
			//printf("mapping term %s\n", ann->goID);
			if((term=getGOTerm(ontology, ann->goID))==HASH_MISS){
				ann++;
				continue;
			}
			if((ann->aspect & aspect)&&
				(evidence & ann->evidence)){

				addReference(term, ontology, markCode);
			}
			ann++;
		}
		markCode++;
		geneName++;
	}
}

char** mapStrings(PyObject* list){
	int len=0;
	char** strings=NULL;
	char** ptr=NULL;
	PyObject* iter=NULL;
	PyObject* item=NULL;
	if(!PyList_Check(list))
		return NULL;
	len=PyList_Size(list)+1;
	strings=malloc(sizeof(char*)*len);
	iter=PyObject_GetIter(list);
	ptr=strings;
	while((item=PyIter_Next(iter))){
		*ptr=PyString_AsString(item);
		Py_DECREF(item);
		ptr++;
	}
	*ptr=NULL;	//sentry
	Py_DECREF(iter);
	return strings;
}

void freeMappedStrings(char** ptr){
	char **tmp=ptr;
	/*while(*ptr){
		free(*ptr);
		ptr++;
	}*/
	free(tmp);
}

double* log_sum_lookup(int size){
	double* data=(double*) malloc(sizeof(double)*(size+1));
	int i=0;
	data[0]=0;
	data[1]=0;
	for(i=2;i<=size;i++)
		data[i]=data[i-1]+log(i);
	return data;
}

double logbin(int n, int r, double* lookup){
	return lookup[n]-lookup[n-r]-lookup[r];

/*	double sum=0;
	int i=0;
	//int t=(n-r+1>2)? n-r+1 : 2;
	for(i=n;i>=n-r+1;--i)
		sum+=log(i);
	for(i=2;i<=r;++i)
		sum-=log(i);
	return sum;*/
}

double binomial(int n, int r, double p, double* lookup){
	if(p==0.0){
		if(r==0.0)
			return 0.0;
		else if(r>0.0)
			return 1.0;
	}else if(p==1.0){
		if(n==r)
			return 1.0;
		else if(n>r)
			return 0.0;
	}
	return exp(logbin(n, r, lookup)+r*log(p)+(n-r)*log(1.0-p));
}
		
double p_value(int nClusterGenes, int nRefGenes, int nGenesMapped, int numRefGenesMapped, double* lookup){
	double sum=0;
	double p=(double)numRefGenesMapped/(double)nRefGenes;
	int i=0;
//	if(p>=1.0 || p<=0.0)
//		printf("p==%f\n",p);
	for(i=nGenesMapped;i<=nClusterGenes;i++)
		sum+=binomial(nClusterGenes,i,p,lookup);
	return sum;
}

PyObject* GOTermFinder(PyObject *self, PyObject* args){
	PyObject* pyGenes=NULL;
	PyObject* pyRefGenes=NULL;
	PyObject* callback=NULL;
	Annotation* annotation=NULL;
	Ontology* ontology=NULL;
	Ontology* slimOntology=NULL;
	int evidence=0;
	int aspect=0;
	int numGenes;
	int numRefGenes=0;
	int progressStep=0;
	double* loglookup=NULL;
	char slimsOnly=0;
	char** genes=NULL;
	char** refGenes=NULL;
	PyObject* result=NULL;
	TermList* list=NULL;
	TermNode* node=NULL;
	//printf("Arg parsing\n");
	if(!PyArg_ParseTuple(args, "OObiiO!O!|OO:GOTermFinder", &pyGenes, &pyRefGenes, &slimsOnly, &evidence, &aspect,
		&go_AnnotationType, &annotation, &go_OntologyType, &ontology,  &slimOntology, &callback)){
		return NULL;
	}
	if(!PyList_Check(pyGenes) || !PyList_Check(pyRefGenes)){
		RAISE_TYPE_ERROR("Need a list");
		return NULL;
	}
	if(callback==Py_None || !PyCallable_Check(callback))
		callback=NULL;
	numGenes=PyList_Size(pyGenes);
	numRefGenes=PyList_Size(pyRefGenes);
	progressStep=(numGenes+numRefGenes)/100;
	//printf("mapping names\n");
	genes=mapStrings(pyGenes);
	refGenes=mapStrings(pyRefGenes);
	prepareGOTerms(ontology);
	//printf("mapping ref\n");
	mapReferenceGenes(refGenes, evidence, aspect, annotation, ontology, callback, MAX(numRefGenes/100,1), 0);
	//printf("finding terms\n");
	list=findTerms(genes, evidence, aspect, slimsOnly, 1, annotation, ontology, slimOntology, callback, MAX(numGenes/100,1), 0);
	node=list->head;
	result=PyDict_New();
	loglookup=log_sum_lookup(numGenes);
	while(node){
		PyObject* geneList=PyList_New(0);
		GOTerm* term=hash_get(ontology->hash, node->id);
		GeneEvidenceNode* ptr=term->mappedGenes.head;
		int numMapped=0;
		double pValue=0;
		while(ptr){
			PyObject* val=Py_BuildValue("s", ptr->geneName);
			PyList_Append(geneList, val);
			Py_DECREF(val);
			numMapped++;
			ptr=ptr->next;
		}
		pValue = p_value(numGenes, numRefGenes, numMapped,  term->numRef, loglookup);
		PyObject* value = Py_BuildValue("Odi", geneList, pValue, term->numRef);
		PyDict_SetItemString(result, term->goID, value);
		Py_DECREF(value);
		Py_DECREF(geneList);
		node=node->next;
	}
	clearTermList(list);
	free(list);
	free(loglookup);
	freeMappedStrings(genes);
	freeMappedStrings(refGenes);
	prepareGOTerms(ontology);
	return result;
}

PyObject* findGOTerms(PyObject* self, PyObject* args){
	PyObject* pyGenes=NULL;
	PyObject* callback=NULL;
	Annotation* annotation=NULL;
	Ontology* ontology=NULL;
	Ontology* slimOntology=NULL;
	int evidence=0;
	int aspect=0;
	char slimsOnly=0;
	char directAnnotation=0;
	char reportEvidences=0;
	char** genes=NULL;
	TermList* list;
	TermNode* node;
	PyObject* result=NULL;
	if(!PyArg_ParseTuple(args, "ObbiibO!O!|OO:findTerms", &pyGenes, &slimsOnly, &directAnnotation, &aspect, &evidence, &reportEvidences,
		&go_AnnotationType, &annotation, &go_OntologyType, &ontology,  &slimOntology, &callback))
		return NULL;
	
	if(!PyList_Check(pyGenes)){
		RAISE_TYPE_ERROR("Need a list");
		return NULL;
	}
	if(callback==Py_None || !PyCallable_Check(callback))
		callback=NULL;
	genes=mapStrings(pyGenes);
	prepareGOTerms(ontology);
	list=findTerms(genes, evidence, aspect, slimsOnly, (char)!directAnnotation, annotation, ontology, slimOntology, callback, MAX(PyList_Size(pyGenes)/100,1),0);
	node=list->head;
	result=PyDict_New();
	while(node){
		PyObject* list=PyList_New(0);
		GOTerm* term=hash_get(ontology->hash, node->id);
		GeneEvidenceNode* ptr=term->mappedGenes.head;
		while(ptr){
			PyObject* val=(reportEvidences)? Py_BuildValue("(si)", ptr->geneName, ptr->evidence) : Py_BuildValue("s", ptr->geneName);
			PyList_Append(list, val);
			Py_DECREF(val);
			ptr=ptr->next;
		}
		PyDict_SetItemString(result, term->goID, list);
		Py_DECREF(list);
		node=node->next;
	}
	clearTermList(list);
	free(list);
	freeMappedStrings(genes);
	return result;
}

PyObject* findGenes(PyObject* self, PyObject* args){
	PyObject* pyTerms=NULL;
	PyObject* result=NULL;
	PyObject* callback=NULL;
	Ontology* ontology=NULL;
	Annotation* annotation=NULL;
	char directAnnotation=0;
	char reportEvidence=0;
	int evidence=0;
	char** terms=NULL;
	char** geneNames=NULL;
	GOTerm* term=NULL;
	AnnRecord* ann=NULL;
	int numGenes=0, i=0;
	char* lastName="No GeNe NaMe 42";
	char** name=NULL;
	TermList* list=NULL;
	TermNode* tnode=NULL;
	GeneEvidenceNode* gnode=NULL;
	hash_table* hash=NULL;

	if(!PyArg_ParseTuple(args, "OibbO!O!|O:findGenes", &pyTerms, &evidence, &reportEvidence, &directAnnotation,
		&go_AnnotationType, &annotation, &go_OntologyType, &ontology, &callback)){
		return NULL;
	}
	if(!PyList_Check(pyTerms)){
		RAISE_TYPE_ERROR("Need a list");
		return NULL;
	}
	if(callback==Py_None || !PyCallable_Check(callback))
		callback=NULL;
	terms=mapStrings(pyTerms);
	prepareGOTerms(ontology);
	ann=annotation->annotation;
	
	while(ann->name){
		numGenes++;
		ann++;
	}
	hash=makeHashTable(PyList_Size(pyTerms));
	name=terms;
	while(*name){
		hash_add(hash, *name, NULL);
		name++;
	}
	ann=annotation->annotation;
	geneNames=malloc(sizeof(char*)*numGenes);
	memset(geneNames, 0, sizeof(char*)*numGenes);
	while(ann->name){
		if(strcmp(lastName, ann->name)!=0){
			lastName=ann->name;
			geneNames[i++]=ann->name;
		}
		ann++;
	}

	list=findTerms(geneNames, evidence, 7, 0, (char)!directAnnotation, annotation, ontology, NULL, callback, numGenes/100, 0);

	tnode=list->head;
	result=PyDict_New();
	while(tnode){
		if((term=getGOTerm(ontology, tnode->id))==HASH_MISS || !hash_has_key(hash, tnode->id)){
			tnode=tnode->next;
			continue;
		}
		gnode=term->mappedGenes.head;
		while(gnode){
			if(PyMapping_HasKeyString(result, gnode->geneName)){
				PyObject* l=PyMapping_GetItemString(result, gnode->geneName);
				PyObject* val=(reportEvidence)? Py_BuildValue("(si)", tnode->id, gnode->evidence) : Py_BuildValue("s", tnode->id);
				PyList_Append(l, val);
				Py_DECREF(l);
				Py_DECREF(val);
			} else{
				PyObject* l=PyList_New(0);
				PyObject* val=(reportEvidence)? Py_BuildValue("(si)", tnode->id, gnode->evidence) : Py_BuildValue("s", tnode->id);
				PyList_Append(l, val);
				PyMapping_SetItemString(result, gnode->geneName, l);
				Py_DECREF(l);
				Py_DECREF(val);
			}
			gnode=gnode->next;
		}
		tnode=tnode->next;
	}
	prepareGOTerms(ontology);
	freeMappedStrings(terms);
	freeHashTable(hash);
	clearTermList(list);
	free(list);
	free(geneNames);
	return result;
}

void collectAllChildrenTerms(GOTerm* term, Ontology* ontology, TermList *list, hash_table* hash){
	TermNode* tnode=NULL;
	if(!hash_has_key(hash, term->goID)){
		addTermNode(list, makeTermNode(term->goID));
		hash_add(hash, term->goID, (void*)42);
	}
	tnode=term->children.head;
	while(tnode){
		term=getGOTerm(ontology, tnode->id);
		if(term!=HASH_MISS)
			collectAllChildrenTerms(term, ontology, list, hash);
		tnode=tnode->next;
	}
}

PyObject* findGenes2(PyObject* self, PyObject* args){
	PyObject* pyTerms=NULL;
	PyObject* result=NULL;
	PyObject* callback=NULL;
	Ontology* ontology=NULL;
	Annotation* annotation=NULL;
	GOTerm* term=NULL;
	char directAnnotation=0;
	char reportEvidence=0;
	int evidence=0;
	char** terms=NULL;
	char** geneNames=NULL;
	hash_table *hash=NULL;
	TermNode *tnode=NULL;
	GeneEvidenceList *gelist=NULL;
	GeneEvidenceNode *genode=NULL;
	AnnRecord *ann=NULL;
	char* lastName="No GeNe NaMe 42";
	char** id=NULL;
	int numTerms=0;
	int numGenes=0;
	TermList *allTerms=makeTermList();
	if(!PyArg_ParseTuple(args, "OibbO!O!|O:findGenes", &pyTerms, &evidence, &reportEvidence, &directAnnotation,
		&go_AnnotationType, &annotation, &go_OntologyType, &ontology, &callback)){
		return NULL;
	}
	if(callback==Py_None || !PyCallable_Check(callback))
		callback=NULL;
	terms=mapStrings(pyTerms);

	id=terms;
	hash=makeHashTable(PyList_Size(pyTerms));
	while(*id){
		hash_add(hash, *id, (void*)42);
		id++;
	}
	if(!directAnnotation){
		hash_table * tmphash=makeHashTable(PyList_Size(pyTerms)*20);
		id=terms;
		while(*id){
			if((term=getGOTerm(ontology, *id))!=HASH_MISS)
				collectAllChildrenTerms(term, ontology, allTerms, tmphash);
			id++;
		}
		tnode=allTerms->head;
		while(tnode){
			numTerms++;
			tnode=tnode->next;
		}
		freeMappedStrings(terms);
		terms=malloc(sizeof(char*)*(numTerms+1));
		tnode=allTerms->head;
		id=terms;
		while(tnode){
			*id=tnode->id;
			id++;
			tnode=tnode->next;
		}
		*id=NULL;
	}
	gelist=makeGeneEvidenceList();
	ann=annotation->annotation;
	while(ann->name){
		if(hash_has_key(hash, ann->goID) && strcmp(ann->name, lastName)){
			addGeneEvidenceNode(gelist, makeGeneEvidenceNode(ann->name, 4095));
			lastName=ann->name;
		}
		ann++;
	}
	genode=gelist->head;
	while(genode){
		numGenes++;
		genode=genode->next;
	}
	geneNames=malloc(sizeof(char*)*(numGenes+1));
	genode=gelist->head;
	id=geneNames;
	while(genode){
		*id=genode->geneName;
		id++;
		genode=genode->next;
	}
	*id=NULL;
	clearTermList(allTerms);
	allTerms=findTerms(geneNames, 4095, 7, 0, (char)!directAnnotation, annotation, ontology, NULL, callback, numGenes/100, 0);
	
	tnode=allTerms->head;
	result=PyDict_New();
	while(tnode){
		if((term=getGOTerm(ontology, tnode->id))==HASH_MISS || !hash_has_key(hash, tnode->id)){
			tnode=tnode->next;
			continue;
		}
		genode=term->mappedGenes.head;
		while(genode){
			if(PyMapping_HasKeyString(result, genode->geneName)){
				PyObject* l=PyMapping_GetItemString(result, genode->geneName);
				PyObject* val=(reportEvidence)? Py_BuildValue("(si)", tnode->id, genode->evidence) : Py_BuildValue("s", tnode->id);
				PyList_Append(l, val);
				Py_DECREF(l);
				Py_DECREF(val);
			} else{
				PyObject* l=PyList_New(0);
				PyObject* val=(reportEvidence)? Py_BuildValue("(si)", tnode->id, genode->evidence) : Py_BuildValue("s", tnode->id);
				PyList_Append(l, val);
				PyMapping_SetItemString(result, genode->geneName, l);
				Py_DECREF(l);
				Py_DECREF(val);
			}
			genode=genode->next;
		}
		tnode=tnode->next;
	}
	prepareGOTerms(ontology);
	freeMappedStrings(terms);
	free(geneNames);
	freeHashTable(hash);
	clearTermList(allTerms);
	clearGeneEvidenceList(gelist);
	free(allTerms);
	free(gelist);
	return result;
}

PyObject* mapToSlims(PyObject* self, PyObject* arg){
	Ontology* ontology=NULL;
	Ontology* slimOntology=NULL;
	char* goId=NULL;
	GOTerm* term=NULL;
	TermList list;
	TermNode* node=NULL;
	PyObject* result=PyList_New(0);
	if(!PyArg_ParseTuple(arg, "sO!O!:mapToSlims", &goId, &go_OntologyType, &ontology, &go_OntologyType, &slimOntology)){
		return NULL;
	}
	prepareGOTerms(ontology);
	memset(&list, 0, sizeof(TermList));
	if((term=getGOTerm(ontology, goId))==HASH_MISS){
		return NULL;
	}
	//printf("mapping slim: %s\n", goId);
	slimMapping(term, ontology, slimOntology, &list);
	//printf("mapping complete\n");
	node=list.head;
	while(node){
		PyObject* val=NULL;
		term=hash_get(ontology->hash, node->id);
		if(term->visited==INDIRECT_MAPPING){
			node=node->next;
			continue;
		}
		val=Py_BuildValue("s",term->goID);
		PyList_Append(result, val);
		Py_DECREF(val);
		node=node->next;
	}
	clearTermList(&list);
	return result;
}

static PyMethodDef go_methods[] = {
	{"parseAnnotation", (PyCFunction)parseAnnotation, METH_VARARGS, "Builds and returns an annotation structure. The argument is a single list containig tuples (geneName, GOId, evidence, aspect) "},
	{"parseGOTerms", (PyCFunction)parseGOTerms, METH_VARARGS, "Builds and returns an ontology structure. The argument is a single list containing tuples (GOid, parents), where parents is a list containg parent GOid's"},
	{"GOTermFinder", (PyCFunction)GOTermFinder, METH_VARARGS, "Maps a list of gene names to relevant terms. Arguments(list of gene names, list of refernece genes, slims only, evidence, aspect, annotation, ontology, slim ontology, callback"},
	{"findTerms", (PyCFunction)findGOTerms, METH_VARARGS, "Maps a list of gene names to terms. Arguments(list of gene names, slims only, direct annotation, evidence, report evidence,  annotation, ontology, slim ontology"},
	{"mapToSlims", (PyCFunction)mapToSlims, METH_VARARGS, "Maps a term to slim terms. Arguments: (term id, full ontology, slim ontology)"},
	{"findGenes", (PyCFunction)findGenes2, METH_VARARGS, "Given a list of GO term id's finds all genes that map to this terms. Arguments: (list of GO term id's, evidence, report evidence, direct annotation only, anotation, ontology "},
    {NULL}  /* Sentinel */
};

#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC
init_GOLib(void) 
{
    PyObject* m;

    go_OntologyType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&go_OntologyType) < 0){
		printf("go type not ready");
        return;
	}

	go_AnnotationType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&go_AnnotationType) < 0){
		printf("annotation type not ready");
        return;
	}

    m = Py_InitModule3("_GOLib", go_methods,
                       "Genetic ontology browser");

    Py_INCREF(&go_OntologyType);
	Py_INCREF(&go_AnnotationType);
    PyModule_AddObject(m, "Ontology", (PyObject *)&go_OntologyType);
	PyModule_AddObject(m, "Annotation", (PyObject *)&go_AnnotationType);
}

