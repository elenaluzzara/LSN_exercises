��

��
8
Const
output"dtype"
valuetensor"
dtypetype

NoOp
C
Placeholder
output"dtype"
dtypetype"
shapeshape:
@
ReadVariableOp
resource
value"dtype"
dtypetype�
�
StatefulPartitionedCall
args2Tin
output2Tout"
Tin
list(type)("
Tout
list(type)("	
ffunc"
configstring "
config_protostring "
executor_typestring �
q
VarHandleOp
resource"
	containerstring "
shared_namestring "
dtypetype"
shapeshape�"serve*2.2.02v2.2.0-rc4-8-g2b96f3662b8��
|
dense_874/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:(*!
shared_namedense_874/kernel
u
$dense_874/kernel/Read/ReadVariableOpReadVariableOpdense_874/kernel*
_output_shapes

:(*
dtype0
t
dense_874/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:(*
shared_namedense_874/bias
m
"dense_874/bias/Read/ReadVariableOpReadVariableOpdense_874/bias*
_output_shapes
:(*
dtype0
|
dense_875/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:((*!
shared_namedense_875/kernel
u
$dense_875/kernel/Read/ReadVariableOpReadVariableOpdense_875/kernel*
_output_shapes

:((*
dtype0
t
dense_875/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:(*
shared_namedense_875/bias
m
"dense_875/bias/Read/ReadVariableOpReadVariableOpdense_875/bias*
_output_shapes
:(*
dtype0
|
dense_876/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:((*!
shared_namedense_876/kernel
u
$dense_876/kernel/Read/ReadVariableOpReadVariableOpdense_876/kernel*
_output_shapes

:((*
dtype0
t
dense_876/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:(*
shared_namedense_876/bias
m
"dense_876/bias/Read/ReadVariableOpReadVariableOpdense_876/bias*
_output_shapes
:(*
dtype0
|
dense_877/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:((*!
shared_namedense_877/kernel
u
$dense_877/kernel/Read/ReadVariableOpReadVariableOpdense_877/kernel*
_output_shapes

:((*
dtype0
t
dense_877/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:(*
shared_namedense_877/bias
m
"dense_877/bias/Read/ReadVariableOpReadVariableOpdense_877/bias*
_output_shapes
:(*
dtype0
|
dense_878/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:((*!
shared_namedense_878/kernel
u
$dense_878/kernel/Read/ReadVariableOpReadVariableOpdense_878/kernel*
_output_shapes

:((*
dtype0
t
dense_878/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:(*
shared_namedense_878/bias
m
"dense_878/bias/Read/ReadVariableOpReadVariableOpdense_878/bias*
_output_shapes
:(*
dtype0
|
dense_879/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:((*!
shared_namedense_879/kernel
u
$dense_879/kernel/Read/ReadVariableOpReadVariableOpdense_879/kernel*
_output_shapes

:((*
dtype0
t
dense_879/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:(*
shared_namedense_879/bias
m
"dense_879/bias/Read/ReadVariableOpReadVariableOpdense_879/bias*
_output_shapes
:(*
dtype0
|
dense_880/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:(*!
shared_namedense_880/kernel
u
$dense_880/kernel/Read/ReadVariableOpReadVariableOpdense_880/kernel*
_output_shapes

:(*
dtype0
t
dense_880/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*
shared_namedense_880/bias
m
"dense_880/bias/Read/ReadVariableOpReadVariableOpdense_880/bias*
_output_shapes
:*
dtype0
n
Adadelta/iterVarHandleOp*
_output_shapes
: *
dtype0	*
shape: *
shared_nameAdadelta/iter
g
!Adadelta/iter/Read/ReadVariableOpReadVariableOpAdadelta/iter*
_output_shapes
: *
dtype0	
p
Adadelta/decayVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nameAdadelta/decay
i
"Adadelta/decay/Read/ReadVariableOpReadVariableOpAdadelta/decay*
_output_shapes
: *
dtype0
�
Adadelta/learning_rateVarHandleOp*
_output_shapes
: *
dtype0*
shape: *'
shared_nameAdadelta/learning_rate
y
*Adadelta/learning_rate/Read/ReadVariableOpReadVariableOpAdadelta/learning_rate*
_output_shapes
: *
dtype0
l
Adadelta/rhoVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nameAdadelta/rho
e
 Adadelta/rho/Read/ReadVariableOpReadVariableOpAdadelta/rho*
_output_shapes
: *
dtype0
^
totalVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nametotal
W
total/Read/ReadVariableOpReadVariableOptotal*
_output_shapes
: *
dtype0
^
countVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_namecount
W
count/Read/ReadVariableOpReadVariableOpcount*
_output_shapes
: *
dtype0
b
total_1VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	total_1
[
total_1/Read/ReadVariableOpReadVariableOptotal_1*
_output_shapes
: *
dtype0
b
count_1VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	count_1
[
count_1/Read/ReadVariableOpReadVariableOpcount_1*
_output_shapes
: *
dtype0
�
$Adadelta/dense_874/kernel/accum_gradVarHandleOp*
_output_shapes
: *
dtype0*
shape
:(*5
shared_name&$Adadelta/dense_874/kernel/accum_grad
�
8Adadelta/dense_874/kernel/accum_grad/Read/ReadVariableOpReadVariableOp$Adadelta/dense_874/kernel/accum_grad*
_output_shapes

:(*
dtype0
�
"Adadelta/dense_874/bias/accum_gradVarHandleOp*
_output_shapes
: *
dtype0*
shape:(*3
shared_name$"Adadelta/dense_874/bias/accum_grad
�
6Adadelta/dense_874/bias/accum_grad/Read/ReadVariableOpReadVariableOp"Adadelta/dense_874/bias/accum_grad*
_output_shapes
:(*
dtype0
�
$Adadelta/dense_875/kernel/accum_gradVarHandleOp*
_output_shapes
: *
dtype0*
shape
:((*5
shared_name&$Adadelta/dense_875/kernel/accum_grad
�
8Adadelta/dense_875/kernel/accum_grad/Read/ReadVariableOpReadVariableOp$Adadelta/dense_875/kernel/accum_grad*
_output_shapes

:((*
dtype0
�
"Adadelta/dense_875/bias/accum_gradVarHandleOp*
_output_shapes
: *
dtype0*
shape:(*3
shared_name$"Adadelta/dense_875/bias/accum_grad
�
6Adadelta/dense_875/bias/accum_grad/Read/ReadVariableOpReadVariableOp"Adadelta/dense_875/bias/accum_grad*
_output_shapes
:(*
dtype0
�
$Adadelta/dense_876/kernel/accum_gradVarHandleOp*
_output_shapes
: *
dtype0*
shape
:((*5
shared_name&$Adadelta/dense_876/kernel/accum_grad
�
8Adadelta/dense_876/kernel/accum_grad/Read/ReadVariableOpReadVariableOp$Adadelta/dense_876/kernel/accum_grad*
_output_shapes

:((*
dtype0
�
"Adadelta/dense_876/bias/accum_gradVarHandleOp*
_output_shapes
: *
dtype0*
shape:(*3
shared_name$"Adadelta/dense_876/bias/accum_grad
�
6Adadelta/dense_876/bias/accum_grad/Read/ReadVariableOpReadVariableOp"Adadelta/dense_876/bias/accum_grad*
_output_shapes
:(*
dtype0
�
$Adadelta/dense_877/kernel/accum_gradVarHandleOp*
_output_shapes
: *
dtype0*
shape
:((*5
shared_name&$Adadelta/dense_877/kernel/accum_grad
�
8Adadelta/dense_877/kernel/accum_grad/Read/ReadVariableOpReadVariableOp$Adadelta/dense_877/kernel/accum_grad*
_output_shapes

:((*
dtype0
�
"Adadelta/dense_877/bias/accum_gradVarHandleOp*
_output_shapes
: *
dtype0*
shape:(*3
shared_name$"Adadelta/dense_877/bias/accum_grad
�
6Adadelta/dense_877/bias/accum_grad/Read/ReadVariableOpReadVariableOp"Adadelta/dense_877/bias/accum_grad*
_output_shapes
:(*
dtype0
�
$Adadelta/dense_878/kernel/accum_gradVarHandleOp*
_output_shapes
: *
dtype0*
shape
:((*5
shared_name&$Adadelta/dense_878/kernel/accum_grad
�
8Adadelta/dense_878/kernel/accum_grad/Read/ReadVariableOpReadVariableOp$Adadelta/dense_878/kernel/accum_grad*
_output_shapes

:((*
dtype0
�
"Adadelta/dense_878/bias/accum_gradVarHandleOp*
_output_shapes
: *
dtype0*
shape:(*3
shared_name$"Adadelta/dense_878/bias/accum_grad
�
6Adadelta/dense_878/bias/accum_grad/Read/ReadVariableOpReadVariableOp"Adadelta/dense_878/bias/accum_grad*
_output_shapes
:(*
dtype0
�
$Adadelta/dense_879/kernel/accum_gradVarHandleOp*
_output_shapes
: *
dtype0*
shape
:((*5
shared_name&$Adadelta/dense_879/kernel/accum_grad
�
8Adadelta/dense_879/kernel/accum_grad/Read/ReadVariableOpReadVariableOp$Adadelta/dense_879/kernel/accum_grad*
_output_shapes

:((*
dtype0
�
"Adadelta/dense_879/bias/accum_gradVarHandleOp*
_output_shapes
: *
dtype0*
shape:(*3
shared_name$"Adadelta/dense_879/bias/accum_grad
�
6Adadelta/dense_879/bias/accum_grad/Read/ReadVariableOpReadVariableOp"Adadelta/dense_879/bias/accum_grad*
_output_shapes
:(*
dtype0
�
$Adadelta/dense_880/kernel/accum_gradVarHandleOp*
_output_shapes
: *
dtype0*
shape
:(*5
shared_name&$Adadelta/dense_880/kernel/accum_grad
�
8Adadelta/dense_880/kernel/accum_grad/Read/ReadVariableOpReadVariableOp$Adadelta/dense_880/kernel/accum_grad*
_output_shapes

:(*
dtype0
�
"Adadelta/dense_880/bias/accum_gradVarHandleOp*
_output_shapes
: *
dtype0*
shape:*3
shared_name$"Adadelta/dense_880/bias/accum_grad
�
6Adadelta/dense_880/bias/accum_grad/Read/ReadVariableOpReadVariableOp"Adadelta/dense_880/bias/accum_grad*
_output_shapes
:*
dtype0
�
#Adadelta/dense_874/kernel/accum_varVarHandleOp*
_output_shapes
: *
dtype0*
shape
:(*4
shared_name%#Adadelta/dense_874/kernel/accum_var
�
7Adadelta/dense_874/kernel/accum_var/Read/ReadVariableOpReadVariableOp#Adadelta/dense_874/kernel/accum_var*
_output_shapes

:(*
dtype0
�
!Adadelta/dense_874/bias/accum_varVarHandleOp*
_output_shapes
: *
dtype0*
shape:(*2
shared_name#!Adadelta/dense_874/bias/accum_var
�
5Adadelta/dense_874/bias/accum_var/Read/ReadVariableOpReadVariableOp!Adadelta/dense_874/bias/accum_var*
_output_shapes
:(*
dtype0
�
#Adadelta/dense_875/kernel/accum_varVarHandleOp*
_output_shapes
: *
dtype0*
shape
:((*4
shared_name%#Adadelta/dense_875/kernel/accum_var
�
7Adadelta/dense_875/kernel/accum_var/Read/ReadVariableOpReadVariableOp#Adadelta/dense_875/kernel/accum_var*
_output_shapes

:((*
dtype0
�
!Adadelta/dense_875/bias/accum_varVarHandleOp*
_output_shapes
: *
dtype0*
shape:(*2
shared_name#!Adadelta/dense_875/bias/accum_var
�
5Adadelta/dense_875/bias/accum_var/Read/ReadVariableOpReadVariableOp!Adadelta/dense_875/bias/accum_var*
_output_shapes
:(*
dtype0
�
#Adadelta/dense_876/kernel/accum_varVarHandleOp*
_output_shapes
: *
dtype0*
shape
:((*4
shared_name%#Adadelta/dense_876/kernel/accum_var
�
7Adadelta/dense_876/kernel/accum_var/Read/ReadVariableOpReadVariableOp#Adadelta/dense_876/kernel/accum_var*
_output_shapes

:((*
dtype0
�
!Adadelta/dense_876/bias/accum_varVarHandleOp*
_output_shapes
: *
dtype0*
shape:(*2
shared_name#!Adadelta/dense_876/bias/accum_var
�
5Adadelta/dense_876/bias/accum_var/Read/ReadVariableOpReadVariableOp!Adadelta/dense_876/bias/accum_var*
_output_shapes
:(*
dtype0
�
#Adadelta/dense_877/kernel/accum_varVarHandleOp*
_output_shapes
: *
dtype0*
shape
:((*4
shared_name%#Adadelta/dense_877/kernel/accum_var
�
7Adadelta/dense_877/kernel/accum_var/Read/ReadVariableOpReadVariableOp#Adadelta/dense_877/kernel/accum_var*
_output_shapes

:((*
dtype0
�
!Adadelta/dense_877/bias/accum_varVarHandleOp*
_output_shapes
: *
dtype0*
shape:(*2
shared_name#!Adadelta/dense_877/bias/accum_var
�
5Adadelta/dense_877/bias/accum_var/Read/ReadVariableOpReadVariableOp!Adadelta/dense_877/bias/accum_var*
_output_shapes
:(*
dtype0
�
#Adadelta/dense_878/kernel/accum_varVarHandleOp*
_output_shapes
: *
dtype0*
shape
:((*4
shared_name%#Adadelta/dense_878/kernel/accum_var
�
7Adadelta/dense_878/kernel/accum_var/Read/ReadVariableOpReadVariableOp#Adadelta/dense_878/kernel/accum_var*
_output_shapes

:((*
dtype0
�
!Adadelta/dense_878/bias/accum_varVarHandleOp*
_output_shapes
: *
dtype0*
shape:(*2
shared_name#!Adadelta/dense_878/bias/accum_var
�
5Adadelta/dense_878/bias/accum_var/Read/ReadVariableOpReadVariableOp!Adadelta/dense_878/bias/accum_var*
_output_shapes
:(*
dtype0
�
#Adadelta/dense_879/kernel/accum_varVarHandleOp*
_output_shapes
: *
dtype0*
shape
:((*4
shared_name%#Adadelta/dense_879/kernel/accum_var
�
7Adadelta/dense_879/kernel/accum_var/Read/ReadVariableOpReadVariableOp#Adadelta/dense_879/kernel/accum_var*
_output_shapes

:((*
dtype0
�
!Adadelta/dense_879/bias/accum_varVarHandleOp*
_output_shapes
: *
dtype0*
shape:(*2
shared_name#!Adadelta/dense_879/bias/accum_var
�
5Adadelta/dense_879/bias/accum_var/Read/ReadVariableOpReadVariableOp!Adadelta/dense_879/bias/accum_var*
_output_shapes
:(*
dtype0
�
#Adadelta/dense_880/kernel/accum_varVarHandleOp*
_output_shapes
: *
dtype0*
shape
:(*4
shared_name%#Adadelta/dense_880/kernel/accum_var
�
7Adadelta/dense_880/kernel/accum_var/Read/ReadVariableOpReadVariableOp#Adadelta/dense_880/kernel/accum_var*
_output_shapes

:(*
dtype0
�
!Adadelta/dense_880/bias/accum_varVarHandleOp*
_output_shapes
: *
dtype0*
shape:*2
shared_name#!Adadelta/dense_880/bias/accum_var
�
5Adadelta/dense_880/bias/accum_var/Read/ReadVariableOpReadVariableOp!Adadelta/dense_880/bias/accum_var*
_output_shapes
:*
dtype0

NoOpNoOp
�O
ConstConst"/device:CPU:0*
_output_shapes
: *
dtype0*�N
value�NB�N B�N
�
layer_with_weights-0
layer-0
layer_with_weights-1
layer-1
layer_with_weights-2
layer-2
layer_with_weights-3
layer-3
layer_with_weights-4
layer-4
layer_with_weights-5
layer-5
layer_with_weights-6
layer-6
	optimizer
	regularization_losses

	variables
trainable_variables
	keras_api

signatures
h

kernel
bias
regularization_losses
	variables
trainable_variables
	keras_api
h

kernel
bias
regularization_losses
	variables
trainable_variables
	keras_api
h

kernel
bias
regularization_losses
	variables
trainable_variables
	keras_api
h

 kernel
!bias
"regularization_losses
#	variables
$trainable_variables
%	keras_api
h

&kernel
'bias
(regularization_losses
)	variables
*trainable_variables
+	keras_api
h

,kernel
-bias
.regularization_losses
/	variables
0trainable_variables
1	keras_api
h

2kernel
3bias
4regularization_losses
5	variables
6trainable_variables
7	keras_api
�
8iter
	9decay
:learning_rate
;rho
accum_grado
accum_gradp
accum_gradq
accum_gradr
accum_grads
accum_gradt 
accum_gradu!
accum_gradv&
accum_gradw'
accum_gradx,
accum_grady-
accum_gradz2
accum_grad{3
accum_grad|	accum_var}	accum_var~	accum_var	accum_var�	accum_var�	accum_var� 	accum_var�!	accum_var�&	accum_var�'	accum_var�,	accum_var�-	accum_var�2	accum_var�3	accum_var�
 
f
0
1
2
3
4
5
 6
!7
&8
'9
,10
-11
212
313
f
0
1
2
3
4
5
 6
!7
&8
'9
,10
-11
212
313
�
<layer_metrics
	regularization_losses
=layer_regularization_losses
>non_trainable_variables

?layers
@metrics

	variables
trainable_variables
 
\Z
VARIABLE_VALUEdense_874/kernel6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUE
XV
VARIABLE_VALUEdense_874/bias4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUE
 

0
1

0
1
�
Alayer_metrics
regularization_losses
Blayer_regularization_losses
Cnon_trainable_variables

Dlayers
Emetrics
	variables
trainable_variables
\Z
VARIABLE_VALUEdense_875/kernel6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUE
XV
VARIABLE_VALUEdense_875/bias4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUE
 

0
1

0
1
�
Flayer_metrics
regularization_losses
Glayer_regularization_losses
Hnon_trainable_variables

Ilayers
Jmetrics
	variables
trainable_variables
\Z
VARIABLE_VALUEdense_876/kernel6layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUE
XV
VARIABLE_VALUEdense_876/bias4layer_with_weights-2/bias/.ATTRIBUTES/VARIABLE_VALUE
 

0
1

0
1
�
Klayer_metrics
regularization_losses
Llayer_regularization_losses
Mnon_trainable_variables

Nlayers
Ometrics
	variables
trainable_variables
\Z
VARIABLE_VALUEdense_877/kernel6layer_with_weights-3/kernel/.ATTRIBUTES/VARIABLE_VALUE
XV
VARIABLE_VALUEdense_877/bias4layer_with_weights-3/bias/.ATTRIBUTES/VARIABLE_VALUE
 

 0
!1

 0
!1
�
Player_metrics
"regularization_losses
Qlayer_regularization_losses
Rnon_trainable_variables

Slayers
Tmetrics
#	variables
$trainable_variables
\Z
VARIABLE_VALUEdense_878/kernel6layer_with_weights-4/kernel/.ATTRIBUTES/VARIABLE_VALUE
XV
VARIABLE_VALUEdense_878/bias4layer_with_weights-4/bias/.ATTRIBUTES/VARIABLE_VALUE
 

&0
'1

&0
'1
�
Ulayer_metrics
(regularization_losses
Vlayer_regularization_losses
Wnon_trainable_variables

Xlayers
Ymetrics
)	variables
*trainable_variables
\Z
VARIABLE_VALUEdense_879/kernel6layer_with_weights-5/kernel/.ATTRIBUTES/VARIABLE_VALUE
XV
VARIABLE_VALUEdense_879/bias4layer_with_weights-5/bias/.ATTRIBUTES/VARIABLE_VALUE
 

,0
-1

,0
-1
�
Zlayer_metrics
.regularization_losses
[layer_regularization_losses
\non_trainable_variables

]layers
^metrics
/	variables
0trainable_variables
\Z
VARIABLE_VALUEdense_880/kernel6layer_with_weights-6/kernel/.ATTRIBUTES/VARIABLE_VALUE
XV
VARIABLE_VALUEdense_880/bias4layer_with_weights-6/bias/.ATTRIBUTES/VARIABLE_VALUE
 

20
31

20
31
�
_layer_metrics
4regularization_losses
`layer_regularization_losses
anon_trainable_variables

blayers
cmetrics
5	variables
6trainable_variables
LJ
VARIABLE_VALUEAdadelta/iter)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUE
NL
VARIABLE_VALUEAdadelta/decay*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUE
^\
VARIABLE_VALUEAdadelta/learning_rate2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUE
JH
VARIABLE_VALUEAdadelta/rho(optimizer/rho/.ATTRIBUTES/VARIABLE_VALUE
 
 
 
1
0
1
2
3
4
5
6

d0
e1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
4
	ftotal
	gcount
h	variables
i	keras_api
D
	jtotal
	kcount
l
_fn_kwargs
m	variables
n	keras_api
OM
VARIABLE_VALUEtotal4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUE
OM
VARIABLE_VALUEcount4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUE

f0
g1

h	variables
QO
VARIABLE_VALUEtotal_14keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUE
QO
VARIABLE_VALUEcount_14keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUE
 

j0
k1

m	variables
��
VARIABLE_VALUE$Adadelta/dense_874/kernel/accum_grad[layer_with_weights-0/kernel/.OPTIMIZER_SLOT/optimizer/accum_grad/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUE"Adadelta/dense_874/bias/accum_gradYlayer_with_weights-0/bias/.OPTIMIZER_SLOT/optimizer/accum_grad/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUE$Adadelta/dense_875/kernel/accum_grad[layer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/accum_grad/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUE"Adadelta/dense_875/bias/accum_gradYlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/accum_grad/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUE$Adadelta/dense_876/kernel/accum_grad[layer_with_weights-2/kernel/.OPTIMIZER_SLOT/optimizer/accum_grad/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUE"Adadelta/dense_876/bias/accum_gradYlayer_with_weights-2/bias/.OPTIMIZER_SLOT/optimizer/accum_grad/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUE$Adadelta/dense_877/kernel/accum_grad[layer_with_weights-3/kernel/.OPTIMIZER_SLOT/optimizer/accum_grad/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUE"Adadelta/dense_877/bias/accum_gradYlayer_with_weights-3/bias/.OPTIMIZER_SLOT/optimizer/accum_grad/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUE$Adadelta/dense_878/kernel/accum_grad[layer_with_weights-4/kernel/.OPTIMIZER_SLOT/optimizer/accum_grad/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUE"Adadelta/dense_878/bias/accum_gradYlayer_with_weights-4/bias/.OPTIMIZER_SLOT/optimizer/accum_grad/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUE$Adadelta/dense_879/kernel/accum_grad[layer_with_weights-5/kernel/.OPTIMIZER_SLOT/optimizer/accum_grad/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUE"Adadelta/dense_879/bias/accum_gradYlayer_with_weights-5/bias/.OPTIMIZER_SLOT/optimizer/accum_grad/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUE$Adadelta/dense_880/kernel/accum_grad[layer_with_weights-6/kernel/.OPTIMIZER_SLOT/optimizer/accum_grad/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUE"Adadelta/dense_880/bias/accum_gradYlayer_with_weights-6/bias/.OPTIMIZER_SLOT/optimizer/accum_grad/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUE#Adadelta/dense_874/kernel/accum_varZlayer_with_weights-0/kernel/.OPTIMIZER_SLOT/optimizer/accum_var/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUE!Adadelta/dense_874/bias/accum_varXlayer_with_weights-0/bias/.OPTIMIZER_SLOT/optimizer/accum_var/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUE#Adadelta/dense_875/kernel/accum_varZlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/accum_var/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUE!Adadelta/dense_875/bias/accum_varXlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/accum_var/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUE#Adadelta/dense_876/kernel/accum_varZlayer_with_weights-2/kernel/.OPTIMIZER_SLOT/optimizer/accum_var/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUE!Adadelta/dense_876/bias/accum_varXlayer_with_weights-2/bias/.OPTIMIZER_SLOT/optimizer/accum_var/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUE#Adadelta/dense_877/kernel/accum_varZlayer_with_weights-3/kernel/.OPTIMIZER_SLOT/optimizer/accum_var/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUE!Adadelta/dense_877/bias/accum_varXlayer_with_weights-3/bias/.OPTIMIZER_SLOT/optimizer/accum_var/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUE#Adadelta/dense_878/kernel/accum_varZlayer_with_weights-4/kernel/.OPTIMIZER_SLOT/optimizer/accum_var/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUE!Adadelta/dense_878/bias/accum_varXlayer_with_weights-4/bias/.OPTIMIZER_SLOT/optimizer/accum_var/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUE#Adadelta/dense_879/kernel/accum_varZlayer_with_weights-5/kernel/.OPTIMIZER_SLOT/optimizer/accum_var/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUE!Adadelta/dense_879/bias/accum_varXlayer_with_weights-5/bias/.OPTIMIZER_SLOT/optimizer/accum_var/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUE#Adadelta/dense_880/kernel/accum_varZlayer_with_weights-6/kernel/.OPTIMIZER_SLOT/optimizer/accum_var/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUE!Adadelta/dense_880/bias/accum_varXlayer_with_weights-6/bias/.OPTIMIZER_SLOT/optimizer/accum_var/.ATTRIBUTES/VARIABLE_VALUE
�
serving_default_dense_874_inputPlaceholder*'
_output_shapes
:���������*
dtype0*
shape:���������
�
StatefulPartitionedCallStatefulPartitionedCallserving_default_dense_874_inputdense_874/kerneldense_874/biasdense_875/kerneldense_875/biasdense_876/kerneldense_876/biasdense_877/kerneldense_877/biasdense_878/kerneldense_878/biasdense_879/kerneldense_879/biasdense_880/kerneldense_880/bias*
Tin
2*
Tout
2*'
_output_shapes
:���������*0
_read_only_resource_inputs
	
**
config_proto

GPU 

CPU2J 8*.
f)R'
%__inference_signature_wrapper_1016964
O
saver_filenamePlaceholder*
_output_shapes
: *
dtype0*
shape: 
�
StatefulPartitionedCall_1StatefulPartitionedCallsaver_filename$dense_874/kernel/Read/ReadVariableOp"dense_874/bias/Read/ReadVariableOp$dense_875/kernel/Read/ReadVariableOp"dense_875/bias/Read/ReadVariableOp$dense_876/kernel/Read/ReadVariableOp"dense_876/bias/Read/ReadVariableOp$dense_877/kernel/Read/ReadVariableOp"dense_877/bias/Read/ReadVariableOp$dense_878/kernel/Read/ReadVariableOp"dense_878/bias/Read/ReadVariableOp$dense_879/kernel/Read/ReadVariableOp"dense_879/bias/Read/ReadVariableOp$dense_880/kernel/Read/ReadVariableOp"dense_880/bias/Read/ReadVariableOp!Adadelta/iter/Read/ReadVariableOp"Adadelta/decay/Read/ReadVariableOp*Adadelta/learning_rate/Read/ReadVariableOp Adadelta/rho/Read/ReadVariableOptotal/Read/ReadVariableOpcount/Read/ReadVariableOptotal_1/Read/ReadVariableOpcount_1/Read/ReadVariableOp8Adadelta/dense_874/kernel/accum_grad/Read/ReadVariableOp6Adadelta/dense_874/bias/accum_grad/Read/ReadVariableOp8Adadelta/dense_875/kernel/accum_grad/Read/ReadVariableOp6Adadelta/dense_875/bias/accum_grad/Read/ReadVariableOp8Adadelta/dense_876/kernel/accum_grad/Read/ReadVariableOp6Adadelta/dense_876/bias/accum_grad/Read/ReadVariableOp8Adadelta/dense_877/kernel/accum_grad/Read/ReadVariableOp6Adadelta/dense_877/bias/accum_grad/Read/ReadVariableOp8Adadelta/dense_878/kernel/accum_grad/Read/ReadVariableOp6Adadelta/dense_878/bias/accum_grad/Read/ReadVariableOp8Adadelta/dense_879/kernel/accum_grad/Read/ReadVariableOp6Adadelta/dense_879/bias/accum_grad/Read/ReadVariableOp8Adadelta/dense_880/kernel/accum_grad/Read/ReadVariableOp6Adadelta/dense_880/bias/accum_grad/Read/ReadVariableOp7Adadelta/dense_874/kernel/accum_var/Read/ReadVariableOp5Adadelta/dense_874/bias/accum_var/Read/ReadVariableOp7Adadelta/dense_875/kernel/accum_var/Read/ReadVariableOp5Adadelta/dense_875/bias/accum_var/Read/ReadVariableOp7Adadelta/dense_876/kernel/accum_var/Read/ReadVariableOp5Adadelta/dense_876/bias/accum_var/Read/ReadVariableOp7Adadelta/dense_877/kernel/accum_var/Read/ReadVariableOp5Adadelta/dense_877/bias/accum_var/Read/ReadVariableOp7Adadelta/dense_878/kernel/accum_var/Read/ReadVariableOp5Adadelta/dense_878/bias/accum_var/Read/ReadVariableOp7Adadelta/dense_879/kernel/accum_var/Read/ReadVariableOp5Adadelta/dense_879/bias/accum_var/Read/ReadVariableOp7Adadelta/dense_880/kernel/accum_var/Read/ReadVariableOp5Adadelta/dense_880/bias/accum_var/Read/ReadVariableOpConst*?
Tin8
624	*
Tout
2*
_output_shapes
: * 
_read_only_resource_inputs
 **
config_proto

GPU 

CPU2J 8*)
f$R"
 __inference__traced_save_1017447
�
StatefulPartitionedCall_2StatefulPartitionedCallsaver_filenamedense_874/kerneldense_874/biasdense_875/kerneldense_875/biasdense_876/kerneldense_876/biasdense_877/kerneldense_877/biasdense_878/kerneldense_878/biasdense_879/kerneldense_879/biasdense_880/kerneldense_880/biasAdadelta/iterAdadelta/decayAdadelta/learning_rateAdadelta/rhototalcounttotal_1count_1$Adadelta/dense_874/kernel/accum_grad"Adadelta/dense_874/bias/accum_grad$Adadelta/dense_875/kernel/accum_grad"Adadelta/dense_875/bias/accum_grad$Adadelta/dense_876/kernel/accum_grad"Adadelta/dense_876/bias/accum_grad$Adadelta/dense_877/kernel/accum_grad"Adadelta/dense_877/bias/accum_grad$Adadelta/dense_878/kernel/accum_grad"Adadelta/dense_878/bias/accum_grad$Adadelta/dense_879/kernel/accum_grad"Adadelta/dense_879/bias/accum_grad$Adadelta/dense_880/kernel/accum_grad"Adadelta/dense_880/bias/accum_grad#Adadelta/dense_874/kernel/accum_var!Adadelta/dense_874/bias/accum_var#Adadelta/dense_875/kernel/accum_var!Adadelta/dense_875/bias/accum_var#Adadelta/dense_876/kernel/accum_var!Adadelta/dense_876/bias/accum_var#Adadelta/dense_877/kernel/accum_var!Adadelta/dense_877/bias/accum_var#Adadelta/dense_878/kernel/accum_var!Adadelta/dense_878/bias/accum_var#Adadelta/dense_879/kernel/accum_var!Adadelta/dense_879/bias/accum_var#Adadelta/dense_880/kernel/accum_var!Adadelta/dense_880/bias/accum_var*>
Tin7
523*
Tout
2*
_output_shapes
: * 
_read_only_resource_inputs
 **
config_proto

GPU 

CPU2J 8*,
f'R%
#__inference__traced_restore_1017609��
�
�
F__inference_dense_877_layer_call_and_return_conditional_losses_1016642

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity��
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:((*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2
MatMul�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:(*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������(2
Reluf
IdentityIdentityRelu:activations:0*
T0*'
_output_shapes
:���������(2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������(:::O K
'
_output_shapes
:���������(
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: 
�
�
F__inference_dense_877_layer_call_and_return_conditional_losses_1017202

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity��
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:((*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2
MatMul�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:(*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������(2
Reluf
IdentityIdentityRelu:activations:0*
T0*'
_output_shapes
:���������(2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������(:::O K
'
_output_shapes
:���������(
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: 
�
�
F__inference_dense_874_layer_call_and_return_conditional_losses_1016561

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity��
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:(*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2
MatMul�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:(*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2	
BiasAddd
IdentityIdentityBiasAdd:output:0*
T0*'
_output_shapes
:���������(2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������:::O K
'
_output_shapes
:���������
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: 
�
�
F__inference_dense_878_layer_call_and_return_conditional_losses_1017222

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity��
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:((*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2
MatMul�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:(*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������(2
Reluf
IdentityIdentityRelu:activations:0*
T0*'
_output_shapes
:���������(2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������(:::O K
'
_output_shapes
:���������(
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: 
�
�
+__inference_dense_879_layer_call_fn_1017251

inputs
unknown
	unknown_0
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*'
_output_shapes
:���������(*$
_read_only_resource_inputs
**
config_proto

GPU 

CPU2J 8*O
fJRH
F__inference_dense_879_layer_call_and_return_conditional_losses_10166962
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������(2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������(::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������(
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: 
�
�
F__inference_dense_880_layer_call_and_return_conditional_losses_1016722

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity��
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:(*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
MatMul�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2	
BiasAddd
IdentityIdentityBiasAdd:output:0*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������(:::O K
'
_output_shapes
:���������(
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: 
�
�
F__inference_dense_875_layer_call_and_return_conditional_losses_1017162

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity��
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:((*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2
MatMul�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:(*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������(2
Reluf
IdentityIdentityRelu:activations:0*
T0*'
_output_shapes
:���������(2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������(:::O K
'
_output_shapes
:���������(
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: 
�
�
F__inference_dense_880_layer_call_and_return_conditional_losses_1017261

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity��
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:(*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
MatMul�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2	
BiasAddd
IdentityIdentityBiasAdd:output:0*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������(:::O K
'
_output_shapes
:���������(
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: 
�5
�
J__inference_sequential_87_layer_call_and_return_conditional_losses_1017015

inputs,
(dense_874_matmul_readvariableop_resource-
)dense_874_biasadd_readvariableop_resource,
(dense_875_matmul_readvariableop_resource-
)dense_875_biasadd_readvariableop_resource,
(dense_876_matmul_readvariableop_resource-
)dense_876_biasadd_readvariableop_resource,
(dense_877_matmul_readvariableop_resource-
)dense_877_biasadd_readvariableop_resource,
(dense_878_matmul_readvariableop_resource-
)dense_878_biasadd_readvariableop_resource,
(dense_879_matmul_readvariableop_resource-
)dense_879_biasadd_readvariableop_resource,
(dense_880_matmul_readvariableop_resource-
)dense_880_biasadd_readvariableop_resource
identity��
dense_874/MatMul/ReadVariableOpReadVariableOp(dense_874_matmul_readvariableop_resource*
_output_shapes

:(*
dtype02!
dense_874/MatMul/ReadVariableOp�
dense_874/MatMulMatMulinputs'dense_874/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2
dense_874/MatMul�
 dense_874/BiasAdd/ReadVariableOpReadVariableOp)dense_874_biasadd_readvariableop_resource*
_output_shapes
:(*
dtype02"
 dense_874/BiasAdd/ReadVariableOp�
dense_874/BiasAddBiasAdddense_874/MatMul:product:0(dense_874/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2
dense_874/BiasAdd�
dense_875/MatMul/ReadVariableOpReadVariableOp(dense_875_matmul_readvariableop_resource*
_output_shapes

:((*
dtype02!
dense_875/MatMul/ReadVariableOp�
dense_875/MatMulMatMuldense_874/BiasAdd:output:0'dense_875/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2
dense_875/MatMul�
 dense_875/BiasAdd/ReadVariableOpReadVariableOp)dense_875_biasadd_readvariableop_resource*
_output_shapes
:(*
dtype02"
 dense_875/BiasAdd/ReadVariableOp�
dense_875/BiasAddBiasAdddense_875/MatMul:product:0(dense_875/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2
dense_875/BiasAddv
dense_875/ReluReludense_875/BiasAdd:output:0*
T0*'
_output_shapes
:���������(2
dense_875/Relu�
dense_876/MatMul/ReadVariableOpReadVariableOp(dense_876_matmul_readvariableop_resource*
_output_shapes

:((*
dtype02!
dense_876/MatMul/ReadVariableOp�
dense_876/MatMulMatMuldense_875/Relu:activations:0'dense_876/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2
dense_876/MatMul�
 dense_876/BiasAdd/ReadVariableOpReadVariableOp)dense_876_biasadd_readvariableop_resource*
_output_shapes
:(*
dtype02"
 dense_876/BiasAdd/ReadVariableOp�
dense_876/BiasAddBiasAdddense_876/MatMul:product:0(dense_876/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2
dense_876/BiasAddv
dense_876/ReluReludense_876/BiasAdd:output:0*
T0*'
_output_shapes
:���������(2
dense_876/Relu�
dense_877/MatMul/ReadVariableOpReadVariableOp(dense_877_matmul_readvariableop_resource*
_output_shapes

:((*
dtype02!
dense_877/MatMul/ReadVariableOp�
dense_877/MatMulMatMuldense_876/Relu:activations:0'dense_877/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2
dense_877/MatMul�
 dense_877/BiasAdd/ReadVariableOpReadVariableOp)dense_877_biasadd_readvariableop_resource*
_output_shapes
:(*
dtype02"
 dense_877/BiasAdd/ReadVariableOp�
dense_877/BiasAddBiasAdddense_877/MatMul:product:0(dense_877/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2
dense_877/BiasAddv
dense_877/ReluReludense_877/BiasAdd:output:0*
T0*'
_output_shapes
:���������(2
dense_877/Relu�
dense_878/MatMul/ReadVariableOpReadVariableOp(dense_878_matmul_readvariableop_resource*
_output_shapes

:((*
dtype02!
dense_878/MatMul/ReadVariableOp�
dense_878/MatMulMatMuldense_877/Relu:activations:0'dense_878/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2
dense_878/MatMul�
 dense_878/BiasAdd/ReadVariableOpReadVariableOp)dense_878_biasadd_readvariableop_resource*
_output_shapes
:(*
dtype02"
 dense_878/BiasAdd/ReadVariableOp�
dense_878/BiasAddBiasAdddense_878/MatMul:product:0(dense_878/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2
dense_878/BiasAddv
dense_878/ReluReludense_878/BiasAdd:output:0*
T0*'
_output_shapes
:���������(2
dense_878/Relu�
dense_879/MatMul/ReadVariableOpReadVariableOp(dense_879_matmul_readvariableop_resource*
_output_shapes

:((*
dtype02!
dense_879/MatMul/ReadVariableOp�
dense_879/MatMulMatMuldense_878/Relu:activations:0'dense_879/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2
dense_879/MatMul�
 dense_879/BiasAdd/ReadVariableOpReadVariableOp)dense_879_biasadd_readvariableop_resource*
_output_shapes
:(*
dtype02"
 dense_879/BiasAdd/ReadVariableOp�
dense_879/BiasAddBiasAdddense_879/MatMul:product:0(dense_879/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2
dense_879/BiasAddv
dense_879/ReluReludense_879/BiasAdd:output:0*
T0*'
_output_shapes
:���������(2
dense_879/Relu�
dense_880/MatMul/ReadVariableOpReadVariableOp(dense_880_matmul_readvariableop_resource*
_output_shapes

:(*
dtype02!
dense_880/MatMul/ReadVariableOp�
dense_880/MatMulMatMuldense_879/Relu:activations:0'dense_880/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
dense_880/MatMul�
 dense_880/BiasAdd/ReadVariableOpReadVariableOp)dense_880_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02"
 dense_880/BiasAdd/ReadVariableOp�
dense_880/BiasAddBiasAdddense_880/MatMul:product:0(dense_880/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
dense_880/BiasAddn
IdentityIdentitydense_880/BiasAdd:output:0*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*^
_input_shapesM
K:���������:::::::::::::::O K
'
_output_shapes
:���������
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :	

_output_shapes
: :


_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: 
��
�
#__inference__traced_restore_1017609
file_prefix%
!assignvariableop_dense_874_kernel%
!assignvariableop_1_dense_874_bias'
#assignvariableop_2_dense_875_kernel%
!assignvariableop_3_dense_875_bias'
#assignvariableop_4_dense_876_kernel%
!assignvariableop_5_dense_876_bias'
#assignvariableop_6_dense_877_kernel%
!assignvariableop_7_dense_877_bias'
#assignvariableop_8_dense_878_kernel%
!assignvariableop_9_dense_878_bias(
$assignvariableop_10_dense_879_kernel&
"assignvariableop_11_dense_879_bias(
$assignvariableop_12_dense_880_kernel&
"assignvariableop_13_dense_880_bias%
!assignvariableop_14_adadelta_iter&
"assignvariableop_15_adadelta_decay.
*assignvariableop_16_adadelta_learning_rate$
 assignvariableop_17_adadelta_rho
assignvariableop_18_total
assignvariableop_19_count
assignvariableop_20_total_1
assignvariableop_21_count_1<
8assignvariableop_22_adadelta_dense_874_kernel_accum_grad:
6assignvariableop_23_adadelta_dense_874_bias_accum_grad<
8assignvariableop_24_adadelta_dense_875_kernel_accum_grad:
6assignvariableop_25_adadelta_dense_875_bias_accum_grad<
8assignvariableop_26_adadelta_dense_876_kernel_accum_grad:
6assignvariableop_27_adadelta_dense_876_bias_accum_grad<
8assignvariableop_28_adadelta_dense_877_kernel_accum_grad:
6assignvariableop_29_adadelta_dense_877_bias_accum_grad<
8assignvariableop_30_adadelta_dense_878_kernel_accum_grad:
6assignvariableop_31_adadelta_dense_878_bias_accum_grad<
8assignvariableop_32_adadelta_dense_879_kernel_accum_grad:
6assignvariableop_33_adadelta_dense_879_bias_accum_grad<
8assignvariableop_34_adadelta_dense_880_kernel_accum_grad:
6assignvariableop_35_adadelta_dense_880_bias_accum_grad;
7assignvariableop_36_adadelta_dense_874_kernel_accum_var9
5assignvariableop_37_adadelta_dense_874_bias_accum_var;
7assignvariableop_38_adadelta_dense_875_kernel_accum_var9
5assignvariableop_39_adadelta_dense_875_bias_accum_var;
7assignvariableop_40_adadelta_dense_876_kernel_accum_var9
5assignvariableop_41_adadelta_dense_876_bias_accum_var;
7assignvariableop_42_adadelta_dense_877_kernel_accum_var9
5assignvariableop_43_adadelta_dense_877_bias_accum_var;
7assignvariableop_44_adadelta_dense_878_kernel_accum_var9
5assignvariableop_45_adadelta_dense_878_bias_accum_var;
7assignvariableop_46_adadelta_dense_879_kernel_accum_var9
5assignvariableop_47_adadelta_dense_879_bias_accum_var;
7assignvariableop_48_adadelta_dense_880_kernel_accum_var9
5assignvariableop_49_adadelta_dense_880_bias_accum_var
identity_51��AssignVariableOp�AssignVariableOp_1�AssignVariableOp_10�AssignVariableOp_11�AssignVariableOp_12�AssignVariableOp_13�AssignVariableOp_14�AssignVariableOp_15�AssignVariableOp_16�AssignVariableOp_17�AssignVariableOp_18�AssignVariableOp_19�AssignVariableOp_2�AssignVariableOp_20�AssignVariableOp_21�AssignVariableOp_22�AssignVariableOp_23�AssignVariableOp_24�AssignVariableOp_25�AssignVariableOp_26�AssignVariableOp_27�AssignVariableOp_28�AssignVariableOp_29�AssignVariableOp_3�AssignVariableOp_30�AssignVariableOp_31�AssignVariableOp_32�AssignVariableOp_33�AssignVariableOp_34�AssignVariableOp_35�AssignVariableOp_36�AssignVariableOp_37�AssignVariableOp_38�AssignVariableOp_39�AssignVariableOp_4�AssignVariableOp_40�AssignVariableOp_41�AssignVariableOp_42�AssignVariableOp_43�AssignVariableOp_44�AssignVariableOp_45�AssignVariableOp_46�AssignVariableOp_47�AssignVariableOp_48�AssignVariableOp_49�AssignVariableOp_5�AssignVariableOp_6�AssignVariableOp_7�AssignVariableOp_8�AssignVariableOp_9�	RestoreV2�RestoreV2_1�
RestoreV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:2*
dtype0*�
value�B�2B6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-2/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-3/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-3/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-4/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-4/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-5/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-5/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-6/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-6/bias/.ATTRIBUTES/VARIABLE_VALUEB)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUEB*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUEB(optimizer/rho/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUEB[layer_with_weights-0/kernel/.OPTIMIZER_SLOT/optimizer/accum_grad/.ATTRIBUTES/VARIABLE_VALUEBYlayer_with_weights-0/bias/.OPTIMIZER_SLOT/optimizer/accum_grad/.ATTRIBUTES/VARIABLE_VALUEB[layer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/accum_grad/.ATTRIBUTES/VARIABLE_VALUEBYlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/accum_grad/.ATTRIBUTES/VARIABLE_VALUEB[layer_with_weights-2/kernel/.OPTIMIZER_SLOT/optimizer/accum_grad/.ATTRIBUTES/VARIABLE_VALUEBYlayer_with_weights-2/bias/.OPTIMIZER_SLOT/optimizer/accum_grad/.ATTRIBUTES/VARIABLE_VALUEB[layer_with_weights-3/kernel/.OPTIMIZER_SLOT/optimizer/accum_grad/.ATTRIBUTES/VARIABLE_VALUEBYlayer_with_weights-3/bias/.OPTIMIZER_SLOT/optimizer/accum_grad/.ATTRIBUTES/VARIABLE_VALUEB[layer_with_weights-4/kernel/.OPTIMIZER_SLOT/optimizer/accum_grad/.ATTRIBUTES/VARIABLE_VALUEBYlayer_with_weights-4/bias/.OPTIMIZER_SLOT/optimizer/accum_grad/.ATTRIBUTES/VARIABLE_VALUEB[layer_with_weights-5/kernel/.OPTIMIZER_SLOT/optimizer/accum_grad/.ATTRIBUTES/VARIABLE_VALUEBYlayer_with_weights-5/bias/.OPTIMIZER_SLOT/optimizer/accum_grad/.ATTRIBUTES/VARIABLE_VALUEB[layer_with_weights-6/kernel/.OPTIMIZER_SLOT/optimizer/accum_grad/.ATTRIBUTES/VARIABLE_VALUEBYlayer_with_weights-6/bias/.OPTIMIZER_SLOT/optimizer/accum_grad/.ATTRIBUTES/VARIABLE_VALUEBZlayer_with_weights-0/kernel/.OPTIMIZER_SLOT/optimizer/accum_var/.ATTRIBUTES/VARIABLE_VALUEBXlayer_with_weights-0/bias/.OPTIMIZER_SLOT/optimizer/accum_var/.ATTRIBUTES/VARIABLE_VALUEBZlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/accum_var/.ATTRIBUTES/VARIABLE_VALUEBXlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/accum_var/.ATTRIBUTES/VARIABLE_VALUEBZlayer_with_weights-2/kernel/.OPTIMIZER_SLOT/optimizer/accum_var/.ATTRIBUTES/VARIABLE_VALUEBXlayer_with_weights-2/bias/.OPTIMIZER_SLOT/optimizer/accum_var/.ATTRIBUTES/VARIABLE_VALUEBZlayer_with_weights-3/kernel/.OPTIMIZER_SLOT/optimizer/accum_var/.ATTRIBUTES/VARIABLE_VALUEBXlayer_with_weights-3/bias/.OPTIMIZER_SLOT/optimizer/accum_var/.ATTRIBUTES/VARIABLE_VALUEBZlayer_with_weights-4/kernel/.OPTIMIZER_SLOT/optimizer/accum_var/.ATTRIBUTES/VARIABLE_VALUEBXlayer_with_weights-4/bias/.OPTIMIZER_SLOT/optimizer/accum_var/.ATTRIBUTES/VARIABLE_VALUEBZlayer_with_weights-5/kernel/.OPTIMIZER_SLOT/optimizer/accum_var/.ATTRIBUTES/VARIABLE_VALUEBXlayer_with_weights-5/bias/.OPTIMIZER_SLOT/optimizer/accum_var/.ATTRIBUTES/VARIABLE_VALUEBZlayer_with_weights-6/kernel/.OPTIMIZER_SLOT/optimizer/accum_var/.ATTRIBUTES/VARIABLE_VALUEBXlayer_with_weights-6/bias/.OPTIMIZER_SLOT/optimizer/accum_var/.ATTRIBUTES/VARIABLE_VALUE2
RestoreV2/tensor_names�
RestoreV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:2*
dtype0*w
valuenBl2B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B 2
RestoreV2/shape_and_slices�
	RestoreV2	RestoreV2file_prefixRestoreV2/tensor_names:output:0#RestoreV2/shape_and_slices:output:0"/device:CPU:0*�
_output_shapes�
�::::::::::::::::::::::::::::::::::::::::::::::::::*@
dtypes6
422	2
	RestoreV2X
IdentityIdentityRestoreV2:tensors:0*
T0*
_output_shapes
:2

Identity�
AssignVariableOpAssignVariableOp!assignvariableop_dense_874_kernelIdentity:output:0*
_output_shapes
 *
dtype02
AssignVariableOp\

Identity_1IdentityRestoreV2:tensors:1*
T0*
_output_shapes
:2

Identity_1�
AssignVariableOp_1AssignVariableOp!assignvariableop_1_dense_874_biasIdentity_1:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_1\

Identity_2IdentityRestoreV2:tensors:2*
T0*
_output_shapes
:2

Identity_2�
AssignVariableOp_2AssignVariableOp#assignvariableop_2_dense_875_kernelIdentity_2:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_2\

Identity_3IdentityRestoreV2:tensors:3*
T0*
_output_shapes
:2

Identity_3�
AssignVariableOp_3AssignVariableOp!assignvariableop_3_dense_875_biasIdentity_3:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_3\

Identity_4IdentityRestoreV2:tensors:4*
T0*
_output_shapes
:2

Identity_4�
AssignVariableOp_4AssignVariableOp#assignvariableop_4_dense_876_kernelIdentity_4:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_4\

Identity_5IdentityRestoreV2:tensors:5*
T0*
_output_shapes
:2

Identity_5�
AssignVariableOp_5AssignVariableOp!assignvariableop_5_dense_876_biasIdentity_5:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_5\

Identity_6IdentityRestoreV2:tensors:6*
T0*
_output_shapes
:2

Identity_6�
AssignVariableOp_6AssignVariableOp#assignvariableop_6_dense_877_kernelIdentity_6:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_6\

Identity_7IdentityRestoreV2:tensors:7*
T0*
_output_shapes
:2

Identity_7�
AssignVariableOp_7AssignVariableOp!assignvariableop_7_dense_877_biasIdentity_7:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_7\

Identity_8IdentityRestoreV2:tensors:8*
T0*
_output_shapes
:2

Identity_8�
AssignVariableOp_8AssignVariableOp#assignvariableop_8_dense_878_kernelIdentity_8:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_8\

Identity_9IdentityRestoreV2:tensors:9*
T0*
_output_shapes
:2

Identity_9�
AssignVariableOp_9AssignVariableOp!assignvariableop_9_dense_878_biasIdentity_9:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_9_
Identity_10IdentityRestoreV2:tensors:10*
T0*
_output_shapes
:2
Identity_10�
AssignVariableOp_10AssignVariableOp$assignvariableop_10_dense_879_kernelIdentity_10:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_10_
Identity_11IdentityRestoreV2:tensors:11*
T0*
_output_shapes
:2
Identity_11�
AssignVariableOp_11AssignVariableOp"assignvariableop_11_dense_879_biasIdentity_11:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_11_
Identity_12IdentityRestoreV2:tensors:12*
T0*
_output_shapes
:2
Identity_12�
AssignVariableOp_12AssignVariableOp$assignvariableop_12_dense_880_kernelIdentity_12:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_12_
Identity_13IdentityRestoreV2:tensors:13*
T0*
_output_shapes
:2
Identity_13�
AssignVariableOp_13AssignVariableOp"assignvariableop_13_dense_880_biasIdentity_13:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_13_
Identity_14IdentityRestoreV2:tensors:14*
T0	*
_output_shapes
:2
Identity_14�
AssignVariableOp_14AssignVariableOp!assignvariableop_14_adadelta_iterIdentity_14:output:0*
_output_shapes
 *
dtype0	2
AssignVariableOp_14_
Identity_15IdentityRestoreV2:tensors:15*
T0*
_output_shapes
:2
Identity_15�
AssignVariableOp_15AssignVariableOp"assignvariableop_15_adadelta_decayIdentity_15:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_15_
Identity_16IdentityRestoreV2:tensors:16*
T0*
_output_shapes
:2
Identity_16�
AssignVariableOp_16AssignVariableOp*assignvariableop_16_adadelta_learning_rateIdentity_16:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_16_
Identity_17IdentityRestoreV2:tensors:17*
T0*
_output_shapes
:2
Identity_17�
AssignVariableOp_17AssignVariableOp assignvariableop_17_adadelta_rhoIdentity_17:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_17_
Identity_18IdentityRestoreV2:tensors:18*
T0*
_output_shapes
:2
Identity_18�
AssignVariableOp_18AssignVariableOpassignvariableop_18_totalIdentity_18:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_18_
Identity_19IdentityRestoreV2:tensors:19*
T0*
_output_shapes
:2
Identity_19�
AssignVariableOp_19AssignVariableOpassignvariableop_19_countIdentity_19:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_19_
Identity_20IdentityRestoreV2:tensors:20*
T0*
_output_shapes
:2
Identity_20�
AssignVariableOp_20AssignVariableOpassignvariableop_20_total_1Identity_20:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_20_
Identity_21IdentityRestoreV2:tensors:21*
T0*
_output_shapes
:2
Identity_21�
AssignVariableOp_21AssignVariableOpassignvariableop_21_count_1Identity_21:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_21_
Identity_22IdentityRestoreV2:tensors:22*
T0*
_output_shapes
:2
Identity_22�
AssignVariableOp_22AssignVariableOp8assignvariableop_22_adadelta_dense_874_kernel_accum_gradIdentity_22:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_22_
Identity_23IdentityRestoreV2:tensors:23*
T0*
_output_shapes
:2
Identity_23�
AssignVariableOp_23AssignVariableOp6assignvariableop_23_adadelta_dense_874_bias_accum_gradIdentity_23:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_23_
Identity_24IdentityRestoreV2:tensors:24*
T0*
_output_shapes
:2
Identity_24�
AssignVariableOp_24AssignVariableOp8assignvariableop_24_adadelta_dense_875_kernel_accum_gradIdentity_24:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_24_
Identity_25IdentityRestoreV2:tensors:25*
T0*
_output_shapes
:2
Identity_25�
AssignVariableOp_25AssignVariableOp6assignvariableop_25_adadelta_dense_875_bias_accum_gradIdentity_25:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_25_
Identity_26IdentityRestoreV2:tensors:26*
T0*
_output_shapes
:2
Identity_26�
AssignVariableOp_26AssignVariableOp8assignvariableop_26_adadelta_dense_876_kernel_accum_gradIdentity_26:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_26_
Identity_27IdentityRestoreV2:tensors:27*
T0*
_output_shapes
:2
Identity_27�
AssignVariableOp_27AssignVariableOp6assignvariableop_27_adadelta_dense_876_bias_accum_gradIdentity_27:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_27_
Identity_28IdentityRestoreV2:tensors:28*
T0*
_output_shapes
:2
Identity_28�
AssignVariableOp_28AssignVariableOp8assignvariableop_28_adadelta_dense_877_kernel_accum_gradIdentity_28:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_28_
Identity_29IdentityRestoreV2:tensors:29*
T0*
_output_shapes
:2
Identity_29�
AssignVariableOp_29AssignVariableOp6assignvariableop_29_adadelta_dense_877_bias_accum_gradIdentity_29:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_29_
Identity_30IdentityRestoreV2:tensors:30*
T0*
_output_shapes
:2
Identity_30�
AssignVariableOp_30AssignVariableOp8assignvariableop_30_adadelta_dense_878_kernel_accum_gradIdentity_30:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_30_
Identity_31IdentityRestoreV2:tensors:31*
T0*
_output_shapes
:2
Identity_31�
AssignVariableOp_31AssignVariableOp6assignvariableop_31_adadelta_dense_878_bias_accum_gradIdentity_31:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_31_
Identity_32IdentityRestoreV2:tensors:32*
T0*
_output_shapes
:2
Identity_32�
AssignVariableOp_32AssignVariableOp8assignvariableop_32_adadelta_dense_879_kernel_accum_gradIdentity_32:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_32_
Identity_33IdentityRestoreV2:tensors:33*
T0*
_output_shapes
:2
Identity_33�
AssignVariableOp_33AssignVariableOp6assignvariableop_33_adadelta_dense_879_bias_accum_gradIdentity_33:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_33_
Identity_34IdentityRestoreV2:tensors:34*
T0*
_output_shapes
:2
Identity_34�
AssignVariableOp_34AssignVariableOp8assignvariableop_34_adadelta_dense_880_kernel_accum_gradIdentity_34:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_34_
Identity_35IdentityRestoreV2:tensors:35*
T0*
_output_shapes
:2
Identity_35�
AssignVariableOp_35AssignVariableOp6assignvariableop_35_adadelta_dense_880_bias_accum_gradIdentity_35:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_35_
Identity_36IdentityRestoreV2:tensors:36*
T0*
_output_shapes
:2
Identity_36�
AssignVariableOp_36AssignVariableOp7assignvariableop_36_adadelta_dense_874_kernel_accum_varIdentity_36:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_36_
Identity_37IdentityRestoreV2:tensors:37*
T0*
_output_shapes
:2
Identity_37�
AssignVariableOp_37AssignVariableOp5assignvariableop_37_adadelta_dense_874_bias_accum_varIdentity_37:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_37_
Identity_38IdentityRestoreV2:tensors:38*
T0*
_output_shapes
:2
Identity_38�
AssignVariableOp_38AssignVariableOp7assignvariableop_38_adadelta_dense_875_kernel_accum_varIdentity_38:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_38_
Identity_39IdentityRestoreV2:tensors:39*
T0*
_output_shapes
:2
Identity_39�
AssignVariableOp_39AssignVariableOp5assignvariableop_39_adadelta_dense_875_bias_accum_varIdentity_39:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_39_
Identity_40IdentityRestoreV2:tensors:40*
T0*
_output_shapes
:2
Identity_40�
AssignVariableOp_40AssignVariableOp7assignvariableop_40_adadelta_dense_876_kernel_accum_varIdentity_40:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_40_
Identity_41IdentityRestoreV2:tensors:41*
T0*
_output_shapes
:2
Identity_41�
AssignVariableOp_41AssignVariableOp5assignvariableop_41_adadelta_dense_876_bias_accum_varIdentity_41:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_41_
Identity_42IdentityRestoreV2:tensors:42*
T0*
_output_shapes
:2
Identity_42�
AssignVariableOp_42AssignVariableOp7assignvariableop_42_adadelta_dense_877_kernel_accum_varIdentity_42:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_42_
Identity_43IdentityRestoreV2:tensors:43*
T0*
_output_shapes
:2
Identity_43�
AssignVariableOp_43AssignVariableOp5assignvariableop_43_adadelta_dense_877_bias_accum_varIdentity_43:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_43_
Identity_44IdentityRestoreV2:tensors:44*
T0*
_output_shapes
:2
Identity_44�
AssignVariableOp_44AssignVariableOp7assignvariableop_44_adadelta_dense_878_kernel_accum_varIdentity_44:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_44_
Identity_45IdentityRestoreV2:tensors:45*
T0*
_output_shapes
:2
Identity_45�
AssignVariableOp_45AssignVariableOp5assignvariableop_45_adadelta_dense_878_bias_accum_varIdentity_45:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_45_
Identity_46IdentityRestoreV2:tensors:46*
T0*
_output_shapes
:2
Identity_46�
AssignVariableOp_46AssignVariableOp7assignvariableop_46_adadelta_dense_879_kernel_accum_varIdentity_46:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_46_
Identity_47IdentityRestoreV2:tensors:47*
T0*
_output_shapes
:2
Identity_47�
AssignVariableOp_47AssignVariableOp5assignvariableop_47_adadelta_dense_879_bias_accum_varIdentity_47:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_47_
Identity_48IdentityRestoreV2:tensors:48*
T0*
_output_shapes
:2
Identity_48�
AssignVariableOp_48AssignVariableOp7assignvariableop_48_adadelta_dense_880_kernel_accum_varIdentity_48:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_48_
Identity_49IdentityRestoreV2:tensors:49*
T0*
_output_shapes
:2
Identity_49�
AssignVariableOp_49AssignVariableOp5assignvariableop_49_adadelta_dense_880_bias_accum_varIdentity_49:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_49�
RestoreV2_1/tensor_namesConst"/device:CPU:0*
_output_shapes
:*
dtype0*1
value(B&B_CHECKPOINTABLE_OBJECT_GRAPH2
RestoreV2_1/tensor_names�
RestoreV2_1/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:*
dtype0*
valueB
B 2
RestoreV2_1/shape_and_slices�
RestoreV2_1	RestoreV2file_prefix!RestoreV2_1/tensor_names:output:0%RestoreV2_1/shape_and_slices:output:0
^RestoreV2"/device:CPU:0*
_output_shapes
:*
dtypes
22
RestoreV2_19
NoOpNoOp"/device:CPU:0*
_output_shapes
 2
NoOp�	
Identity_50Identityfile_prefix^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_11^AssignVariableOp_12^AssignVariableOp_13^AssignVariableOp_14^AssignVariableOp_15^AssignVariableOp_16^AssignVariableOp_17^AssignVariableOp_18^AssignVariableOp_19^AssignVariableOp_2^AssignVariableOp_20^AssignVariableOp_21^AssignVariableOp_22^AssignVariableOp_23^AssignVariableOp_24^AssignVariableOp_25^AssignVariableOp_26^AssignVariableOp_27^AssignVariableOp_28^AssignVariableOp_29^AssignVariableOp_3^AssignVariableOp_30^AssignVariableOp_31^AssignVariableOp_32^AssignVariableOp_33^AssignVariableOp_34^AssignVariableOp_35^AssignVariableOp_36^AssignVariableOp_37^AssignVariableOp_38^AssignVariableOp_39^AssignVariableOp_4^AssignVariableOp_40^AssignVariableOp_41^AssignVariableOp_42^AssignVariableOp_43^AssignVariableOp_44^AssignVariableOp_45^AssignVariableOp_46^AssignVariableOp_47^AssignVariableOp_48^AssignVariableOp_49^AssignVariableOp_5^AssignVariableOp_6^AssignVariableOp_7^AssignVariableOp_8^AssignVariableOp_9^NoOp"/device:CPU:0*
T0*
_output_shapes
: 2
Identity_50�	
Identity_51IdentityIdentity_50:output:0^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_11^AssignVariableOp_12^AssignVariableOp_13^AssignVariableOp_14^AssignVariableOp_15^AssignVariableOp_16^AssignVariableOp_17^AssignVariableOp_18^AssignVariableOp_19^AssignVariableOp_2^AssignVariableOp_20^AssignVariableOp_21^AssignVariableOp_22^AssignVariableOp_23^AssignVariableOp_24^AssignVariableOp_25^AssignVariableOp_26^AssignVariableOp_27^AssignVariableOp_28^AssignVariableOp_29^AssignVariableOp_3^AssignVariableOp_30^AssignVariableOp_31^AssignVariableOp_32^AssignVariableOp_33^AssignVariableOp_34^AssignVariableOp_35^AssignVariableOp_36^AssignVariableOp_37^AssignVariableOp_38^AssignVariableOp_39^AssignVariableOp_4^AssignVariableOp_40^AssignVariableOp_41^AssignVariableOp_42^AssignVariableOp_43^AssignVariableOp_44^AssignVariableOp_45^AssignVariableOp_46^AssignVariableOp_47^AssignVariableOp_48^AssignVariableOp_49^AssignVariableOp_5^AssignVariableOp_6^AssignVariableOp_7^AssignVariableOp_8^AssignVariableOp_9
^RestoreV2^RestoreV2_1*
T0*
_output_shapes
: 2
Identity_51"#
identity_51Identity_51:output:0*�
_input_shapes�
�: ::::::::::::::::::::::::::::::::::::::::::::::::::2$
AssignVariableOpAssignVariableOp2(
AssignVariableOp_1AssignVariableOp_12*
AssignVariableOp_10AssignVariableOp_102*
AssignVariableOp_11AssignVariableOp_112*
AssignVariableOp_12AssignVariableOp_122*
AssignVariableOp_13AssignVariableOp_132*
AssignVariableOp_14AssignVariableOp_142*
AssignVariableOp_15AssignVariableOp_152*
AssignVariableOp_16AssignVariableOp_162*
AssignVariableOp_17AssignVariableOp_172*
AssignVariableOp_18AssignVariableOp_182*
AssignVariableOp_19AssignVariableOp_192(
AssignVariableOp_2AssignVariableOp_22*
AssignVariableOp_20AssignVariableOp_202*
AssignVariableOp_21AssignVariableOp_212*
AssignVariableOp_22AssignVariableOp_222*
AssignVariableOp_23AssignVariableOp_232*
AssignVariableOp_24AssignVariableOp_242*
AssignVariableOp_25AssignVariableOp_252*
AssignVariableOp_26AssignVariableOp_262*
AssignVariableOp_27AssignVariableOp_272*
AssignVariableOp_28AssignVariableOp_282*
AssignVariableOp_29AssignVariableOp_292(
AssignVariableOp_3AssignVariableOp_32*
AssignVariableOp_30AssignVariableOp_302*
AssignVariableOp_31AssignVariableOp_312*
AssignVariableOp_32AssignVariableOp_322*
AssignVariableOp_33AssignVariableOp_332*
AssignVariableOp_34AssignVariableOp_342*
AssignVariableOp_35AssignVariableOp_352*
AssignVariableOp_36AssignVariableOp_362*
AssignVariableOp_37AssignVariableOp_372*
AssignVariableOp_38AssignVariableOp_382*
AssignVariableOp_39AssignVariableOp_392(
AssignVariableOp_4AssignVariableOp_42*
AssignVariableOp_40AssignVariableOp_402*
AssignVariableOp_41AssignVariableOp_412*
AssignVariableOp_42AssignVariableOp_422*
AssignVariableOp_43AssignVariableOp_432*
AssignVariableOp_44AssignVariableOp_442*
AssignVariableOp_45AssignVariableOp_452*
AssignVariableOp_46AssignVariableOp_462*
AssignVariableOp_47AssignVariableOp_472*
AssignVariableOp_48AssignVariableOp_482*
AssignVariableOp_49AssignVariableOp_492(
AssignVariableOp_5AssignVariableOp_52(
AssignVariableOp_6AssignVariableOp_62(
AssignVariableOp_7AssignVariableOp_72(
AssignVariableOp_8AssignVariableOp_82(
AssignVariableOp_9AssignVariableOp_92
	RestoreV2	RestoreV22
RestoreV2_1RestoreV2_1:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix:

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :	

_output_shapes
: :


_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: : 

_output_shapes
: :!

_output_shapes
: :"

_output_shapes
: :#

_output_shapes
: :$

_output_shapes
: :%

_output_shapes
: :&

_output_shapes
: :'

_output_shapes
: :(

_output_shapes
: :)

_output_shapes
: :*

_output_shapes
: :+

_output_shapes
: :,

_output_shapes
: :-

_output_shapes
: :.

_output_shapes
: :/

_output_shapes
: :0

_output_shapes
: :1

_output_shapes
: :2

_output_shapes
: 
�
�
+__inference_dense_875_layer_call_fn_1017171

inputs
unknown
	unknown_0
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*'
_output_shapes
:���������(*$
_read_only_resource_inputs
**
config_proto

GPU 

CPU2J 8*O
fJRH
F__inference_dense_875_layer_call_and_return_conditional_losses_10165882
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������(2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������(::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������(
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: 
�5
�
J__inference_sequential_87_layer_call_and_return_conditional_losses_1017066

inputs,
(dense_874_matmul_readvariableop_resource-
)dense_874_biasadd_readvariableop_resource,
(dense_875_matmul_readvariableop_resource-
)dense_875_biasadd_readvariableop_resource,
(dense_876_matmul_readvariableop_resource-
)dense_876_biasadd_readvariableop_resource,
(dense_877_matmul_readvariableop_resource-
)dense_877_biasadd_readvariableop_resource,
(dense_878_matmul_readvariableop_resource-
)dense_878_biasadd_readvariableop_resource,
(dense_879_matmul_readvariableop_resource-
)dense_879_biasadd_readvariableop_resource,
(dense_880_matmul_readvariableop_resource-
)dense_880_biasadd_readvariableop_resource
identity��
dense_874/MatMul/ReadVariableOpReadVariableOp(dense_874_matmul_readvariableop_resource*
_output_shapes

:(*
dtype02!
dense_874/MatMul/ReadVariableOp�
dense_874/MatMulMatMulinputs'dense_874/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2
dense_874/MatMul�
 dense_874/BiasAdd/ReadVariableOpReadVariableOp)dense_874_biasadd_readvariableop_resource*
_output_shapes
:(*
dtype02"
 dense_874/BiasAdd/ReadVariableOp�
dense_874/BiasAddBiasAdddense_874/MatMul:product:0(dense_874/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2
dense_874/BiasAdd�
dense_875/MatMul/ReadVariableOpReadVariableOp(dense_875_matmul_readvariableop_resource*
_output_shapes

:((*
dtype02!
dense_875/MatMul/ReadVariableOp�
dense_875/MatMulMatMuldense_874/BiasAdd:output:0'dense_875/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2
dense_875/MatMul�
 dense_875/BiasAdd/ReadVariableOpReadVariableOp)dense_875_biasadd_readvariableop_resource*
_output_shapes
:(*
dtype02"
 dense_875/BiasAdd/ReadVariableOp�
dense_875/BiasAddBiasAdddense_875/MatMul:product:0(dense_875/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2
dense_875/BiasAddv
dense_875/ReluReludense_875/BiasAdd:output:0*
T0*'
_output_shapes
:���������(2
dense_875/Relu�
dense_876/MatMul/ReadVariableOpReadVariableOp(dense_876_matmul_readvariableop_resource*
_output_shapes

:((*
dtype02!
dense_876/MatMul/ReadVariableOp�
dense_876/MatMulMatMuldense_875/Relu:activations:0'dense_876/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2
dense_876/MatMul�
 dense_876/BiasAdd/ReadVariableOpReadVariableOp)dense_876_biasadd_readvariableop_resource*
_output_shapes
:(*
dtype02"
 dense_876/BiasAdd/ReadVariableOp�
dense_876/BiasAddBiasAdddense_876/MatMul:product:0(dense_876/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2
dense_876/BiasAddv
dense_876/ReluReludense_876/BiasAdd:output:0*
T0*'
_output_shapes
:���������(2
dense_876/Relu�
dense_877/MatMul/ReadVariableOpReadVariableOp(dense_877_matmul_readvariableop_resource*
_output_shapes

:((*
dtype02!
dense_877/MatMul/ReadVariableOp�
dense_877/MatMulMatMuldense_876/Relu:activations:0'dense_877/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2
dense_877/MatMul�
 dense_877/BiasAdd/ReadVariableOpReadVariableOp)dense_877_biasadd_readvariableop_resource*
_output_shapes
:(*
dtype02"
 dense_877/BiasAdd/ReadVariableOp�
dense_877/BiasAddBiasAdddense_877/MatMul:product:0(dense_877/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2
dense_877/BiasAddv
dense_877/ReluReludense_877/BiasAdd:output:0*
T0*'
_output_shapes
:���������(2
dense_877/Relu�
dense_878/MatMul/ReadVariableOpReadVariableOp(dense_878_matmul_readvariableop_resource*
_output_shapes

:((*
dtype02!
dense_878/MatMul/ReadVariableOp�
dense_878/MatMulMatMuldense_877/Relu:activations:0'dense_878/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2
dense_878/MatMul�
 dense_878/BiasAdd/ReadVariableOpReadVariableOp)dense_878_biasadd_readvariableop_resource*
_output_shapes
:(*
dtype02"
 dense_878/BiasAdd/ReadVariableOp�
dense_878/BiasAddBiasAdddense_878/MatMul:product:0(dense_878/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2
dense_878/BiasAddv
dense_878/ReluReludense_878/BiasAdd:output:0*
T0*'
_output_shapes
:���������(2
dense_878/Relu�
dense_879/MatMul/ReadVariableOpReadVariableOp(dense_879_matmul_readvariableop_resource*
_output_shapes

:((*
dtype02!
dense_879/MatMul/ReadVariableOp�
dense_879/MatMulMatMuldense_878/Relu:activations:0'dense_879/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2
dense_879/MatMul�
 dense_879/BiasAdd/ReadVariableOpReadVariableOp)dense_879_biasadd_readvariableop_resource*
_output_shapes
:(*
dtype02"
 dense_879/BiasAdd/ReadVariableOp�
dense_879/BiasAddBiasAdddense_879/MatMul:product:0(dense_879/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2
dense_879/BiasAddv
dense_879/ReluReludense_879/BiasAdd:output:0*
T0*'
_output_shapes
:���������(2
dense_879/Relu�
dense_880/MatMul/ReadVariableOpReadVariableOp(dense_880_matmul_readvariableop_resource*
_output_shapes

:(*
dtype02!
dense_880/MatMul/ReadVariableOp�
dense_880/MatMulMatMuldense_879/Relu:activations:0'dense_880/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
dense_880/MatMul�
 dense_880/BiasAdd/ReadVariableOpReadVariableOp)dense_880_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02"
 dense_880/BiasAdd/ReadVariableOp�
dense_880/BiasAddBiasAdddense_880/MatMul:product:0(dense_880/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
dense_880/BiasAddn
IdentityIdentitydense_880/BiasAdd:output:0*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*^
_input_shapesM
K:���������:::::::::::::::O K
'
_output_shapes
:���������
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :	

_output_shapes
: :


_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: 
�'
�
J__inference_sequential_87_layer_call_and_return_conditional_losses_1016739
dense_874_input
dense_874_1016572
dense_874_1016574
dense_875_1016599
dense_875_1016601
dense_876_1016626
dense_876_1016628
dense_877_1016653
dense_877_1016655
dense_878_1016680
dense_878_1016682
dense_879_1016707
dense_879_1016709
dense_880_1016733
dense_880_1016735
identity��!dense_874/StatefulPartitionedCall�!dense_875/StatefulPartitionedCall�!dense_876/StatefulPartitionedCall�!dense_877/StatefulPartitionedCall�!dense_878/StatefulPartitionedCall�!dense_879/StatefulPartitionedCall�!dense_880/StatefulPartitionedCall�
!dense_874/StatefulPartitionedCallStatefulPartitionedCalldense_874_inputdense_874_1016572dense_874_1016574*
Tin
2*
Tout
2*'
_output_shapes
:���������(*$
_read_only_resource_inputs
**
config_proto

GPU 

CPU2J 8*O
fJRH
F__inference_dense_874_layer_call_and_return_conditional_losses_10165612#
!dense_874/StatefulPartitionedCall�
!dense_875/StatefulPartitionedCallStatefulPartitionedCall*dense_874/StatefulPartitionedCall:output:0dense_875_1016599dense_875_1016601*
Tin
2*
Tout
2*'
_output_shapes
:���������(*$
_read_only_resource_inputs
**
config_proto

GPU 

CPU2J 8*O
fJRH
F__inference_dense_875_layer_call_and_return_conditional_losses_10165882#
!dense_875/StatefulPartitionedCall�
!dense_876/StatefulPartitionedCallStatefulPartitionedCall*dense_875/StatefulPartitionedCall:output:0dense_876_1016626dense_876_1016628*
Tin
2*
Tout
2*'
_output_shapes
:���������(*$
_read_only_resource_inputs
**
config_proto

GPU 

CPU2J 8*O
fJRH
F__inference_dense_876_layer_call_and_return_conditional_losses_10166152#
!dense_876/StatefulPartitionedCall�
!dense_877/StatefulPartitionedCallStatefulPartitionedCall*dense_876/StatefulPartitionedCall:output:0dense_877_1016653dense_877_1016655*
Tin
2*
Tout
2*'
_output_shapes
:���������(*$
_read_only_resource_inputs
**
config_proto

GPU 

CPU2J 8*O
fJRH
F__inference_dense_877_layer_call_and_return_conditional_losses_10166422#
!dense_877/StatefulPartitionedCall�
!dense_878/StatefulPartitionedCallStatefulPartitionedCall*dense_877/StatefulPartitionedCall:output:0dense_878_1016680dense_878_1016682*
Tin
2*
Tout
2*'
_output_shapes
:���������(*$
_read_only_resource_inputs
**
config_proto

GPU 

CPU2J 8*O
fJRH
F__inference_dense_878_layer_call_and_return_conditional_losses_10166692#
!dense_878/StatefulPartitionedCall�
!dense_879/StatefulPartitionedCallStatefulPartitionedCall*dense_878/StatefulPartitionedCall:output:0dense_879_1016707dense_879_1016709*
Tin
2*
Tout
2*'
_output_shapes
:���������(*$
_read_only_resource_inputs
**
config_proto

GPU 

CPU2J 8*O
fJRH
F__inference_dense_879_layer_call_and_return_conditional_losses_10166962#
!dense_879/StatefulPartitionedCall�
!dense_880/StatefulPartitionedCallStatefulPartitionedCall*dense_879/StatefulPartitionedCall:output:0dense_880_1016733dense_880_1016735*
Tin
2*
Tout
2*'
_output_shapes
:���������*$
_read_only_resource_inputs
**
config_proto

GPU 

CPU2J 8*O
fJRH
F__inference_dense_880_layer_call_and_return_conditional_losses_10167222#
!dense_880/StatefulPartitionedCall�
IdentityIdentity*dense_880/StatefulPartitionedCall:output:0"^dense_874/StatefulPartitionedCall"^dense_875/StatefulPartitionedCall"^dense_876/StatefulPartitionedCall"^dense_877/StatefulPartitionedCall"^dense_878/StatefulPartitionedCall"^dense_879/StatefulPartitionedCall"^dense_880/StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*^
_input_shapesM
K:���������::::::::::::::2F
!dense_874/StatefulPartitionedCall!dense_874/StatefulPartitionedCall2F
!dense_875/StatefulPartitionedCall!dense_875/StatefulPartitionedCall2F
!dense_876/StatefulPartitionedCall!dense_876/StatefulPartitionedCall2F
!dense_877/StatefulPartitionedCall!dense_877/StatefulPartitionedCall2F
!dense_878/StatefulPartitionedCall!dense_878/StatefulPartitionedCall2F
!dense_879/StatefulPartitionedCall!dense_879/StatefulPartitionedCall2F
!dense_880/StatefulPartitionedCall!dense_880/StatefulPartitionedCall:X T
'
_output_shapes
:���������
)
_user_specified_namedense_874_input:

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :	

_output_shapes
: :


_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: 
�
�
F__inference_dense_879_layer_call_and_return_conditional_losses_1016696

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity��
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:((*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2
MatMul�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:(*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������(2
Reluf
IdentityIdentityRelu:activations:0*
T0*'
_output_shapes
:���������(2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������(:::O K
'
_output_shapes
:���������(
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: 
�
�
/__inference_sequential_87_layer_call_fn_1016923
dense_874_input
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
	unknown_9

unknown_10

unknown_11

unknown_12
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCalldense_874_inputunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12*
Tin
2*
Tout
2*'
_output_shapes
:���������*0
_read_only_resource_inputs
	
**
config_proto

GPU 

CPU2J 8*S
fNRL
J__inference_sequential_87_layer_call_and_return_conditional_losses_10168922
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*^
_input_shapesM
K:���������::::::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:X T
'
_output_shapes
:���������
)
_user_specified_namedense_874_input:

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :	

_output_shapes
: :


_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: 
�'
�
J__inference_sequential_87_layer_call_and_return_conditional_losses_1016778
dense_874_input
dense_874_1016742
dense_874_1016744
dense_875_1016747
dense_875_1016749
dense_876_1016752
dense_876_1016754
dense_877_1016757
dense_877_1016759
dense_878_1016762
dense_878_1016764
dense_879_1016767
dense_879_1016769
dense_880_1016772
dense_880_1016774
identity��!dense_874/StatefulPartitionedCall�!dense_875/StatefulPartitionedCall�!dense_876/StatefulPartitionedCall�!dense_877/StatefulPartitionedCall�!dense_878/StatefulPartitionedCall�!dense_879/StatefulPartitionedCall�!dense_880/StatefulPartitionedCall�
!dense_874/StatefulPartitionedCallStatefulPartitionedCalldense_874_inputdense_874_1016742dense_874_1016744*
Tin
2*
Tout
2*'
_output_shapes
:���������(*$
_read_only_resource_inputs
**
config_proto

GPU 

CPU2J 8*O
fJRH
F__inference_dense_874_layer_call_and_return_conditional_losses_10165612#
!dense_874/StatefulPartitionedCall�
!dense_875/StatefulPartitionedCallStatefulPartitionedCall*dense_874/StatefulPartitionedCall:output:0dense_875_1016747dense_875_1016749*
Tin
2*
Tout
2*'
_output_shapes
:���������(*$
_read_only_resource_inputs
**
config_proto

GPU 

CPU2J 8*O
fJRH
F__inference_dense_875_layer_call_and_return_conditional_losses_10165882#
!dense_875/StatefulPartitionedCall�
!dense_876/StatefulPartitionedCallStatefulPartitionedCall*dense_875/StatefulPartitionedCall:output:0dense_876_1016752dense_876_1016754*
Tin
2*
Tout
2*'
_output_shapes
:���������(*$
_read_only_resource_inputs
**
config_proto

GPU 

CPU2J 8*O
fJRH
F__inference_dense_876_layer_call_and_return_conditional_losses_10166152#
!dense_876/StatefulPartitionedCall�
!dense_877/StatefulPartitionedCallStatefulPartitionedCall*dense_876/StatefulPartitionedCall:output:0dense_877_1016757dense_877_1016759*
Tin
2*
Tout
2*'
_output_shapes
:���������(*$
_read_only_resource_inputs
**
config_proto

GPU 

CPU2J 8*O
fJRH
F__inference_dense_877_layer_call_and_return_conditional_losses_10166422#
!dense_877/StatefulPartitionedCall�
!dense_878/StatefulPartitionedCallStatefulPartitionedCall*dense_877/StatefulPartitionedCall:output:0dense_878_1016762dense_878_1016764*
Tin
2*
Tout
2*'
_output_shapes
:���������(*$
_read_only_resource_inputs
**
config_proto

GPU 

CPU2J 8*O
fJRH
F__inference_dense_878_layer_call_and_return_conditional_losses_10166692#
!dense_878/StatefulPartitionedCall�
!dense_879/StatefulPartitionedCallStatefulPartitionedCall*dense_878/StatefulPartitionedCall:output:0dense_879_1016767dense_879_1016769*
Tin
2*
Tout
2*'
_output_shapes
:���������(*$
_read_only_resource_inputs
**
config_proto

GPU 

CPU2J 8*O
fJRH
F__inference_dense_879_layer_call_and_return_conditional_losses_10166962#
!dense_879/StatefulPartitionedCall�
!dense_880/StatefulPartitionedCallStatefulPartitionedCall*dense_879/StatefulPartitionedCall:output:0dense_880_1016772dense_880_1016774*
Tin
2*
Tout
2*'
_output_shapes
:���������*$
_read_only_resource_inputs
**
config_proto

GPU 

CPU2J 8*O
fJRH
F__inference_dense_880_layer_call_and_return_conditional_losses_10167222#
!dense_880/StatefulPartitionedCall�
IdentityIdentity*dense_880/StatefulPartitionedCall:output:0"^dense_874/StatefulPartitionedCall"^dense_875/StatefulPartitionedCall"^dense_876/StatefulPartitionedCall"^dense_877/StatefulPartitionedCall"^dense_878/StatefulPartitionedCall"^dense_879/StatefulPartitionedCall"^dense_880/StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*^
_input_shapesM
K:���������::::::::::::::2F
!dense_874/StatefulPartitionedCall!dense_874/StatefulPartitionedCall2F
!dense_875/StatefulPartitionedCall!dense_875/StatefulPartitionedCall2F
!dense_876/StatefulPartitionedCall!dense_876/StatefulPartitionedCall2F
!dense_877/StatefulPartitionedCall!dense_877/StatefulPartitionedCall2F
!dense_878/StatefulPartitionedCall!dense_878/StatefulPartitionedCall2F
!dense_879/StatefulPartitionedCall!dense_879/StatefulPartitionedCall2F
!dense_880/StatefulPartitionedCall!dense_880/StatefulPartitionedCall:X T
'
_output_shapes
:���������
)
_user_specified_namedense_874_input:

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :	

_output_shapes
: :


_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: 
�
�
%__inference_signature_wrapper_1016964
dense_874_input
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
	unknown_9

unknown_10

unknown_11

unknown_12
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCalldense_874_inputunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12*
Tin
2*
Tout
2*'
_output_shapes
:���������*0
_read_only_resource_inputs
	
**
config_proto

GPU 

CPU2J 8*+
f&R$
"__inference__wrapped_model_10165472
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*^
_input_shapesM
K:���������::::::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:X T
'
_output_shapes
:���������
)
_user_specified_namedense_874_input:

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :	

_output_shapes
: :


_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: 
�
�
+__inference_dense_877_layer_call_fn_1017211

inputs
unknown
	unknown_0
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*'
_output_shapes
:���������(*$
_read_only_resource_inputs
**
config_proto

GPU 

CPU2J 8*O
fJRH
F__inference_dense_877_layer_call_and_return_conditional_losses_10166422
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������(2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������(::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������(
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: 
�
�
F__inference_dense_874_layer_call_and_return_conditional_losses_1017142

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity��
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:(*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2
MatMul�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:(*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2	
BiasAddd
IdentityIdentityBiasAdd:output:0*
T0*'
_output_shapes
:���������(2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������:::O K
'
_output_shapes
:���������
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: 
�
�
+__inference_dense_878_layer_call_fn_1017231

inputs
unknown
	unknown_0
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*'
_output_shapes
:���������(*$
_read_only_resource_inputs
**
config_proto

GPU 

CPU2J 8*O
fJRH
F__inference_dense_878_layer_call_and_return_conditional_losses_10166692
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������(2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������(::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������(
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: 
�
�
F__inference_dense_879_layer_call_and_return_conditional_losses_1017242

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity��
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:((*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2
MatMul�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:(*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������(2
Reluf
IdentityIdentityRelu:activations:0*
T0*'
_output_shapes
:���������(2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������(:::O K
'
_output_shapes
:���������(
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: 
�
�
/__inference_sequential_87_layer_call_fn_1017099

inputs
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
	unknown_9

unknown_10

unknown_11

unknown_12
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12*
Tin
2*
Tout
2*'
_output_shapes
:���������*0
_read_only_resource_inputs
	
**
config_proto

GPU 

CPU2J 8*S
fNRL
J__inference_sequential_87_layer_call_and_return_conditional_losses_10168202
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*^
_input_shapesM
K:���������::::::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :	

_output_shapes
: :


_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: 
�
�
+__inference_dense_874_layer_call_fn_1017151

inputs
unknown
	unknown_0
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*'
_output_shapes
:���������(*$
_read_only_resource_inputs
**
config_proto

GPU 

CPU2J 8*O
fJRH
F__inference_dense_874_layer_call_and_return_conditional_losses_10165612
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������(2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: 
�
�
+__inference_dense_880_layer_call_fn_1017270

inputs
unknown
	unknown_0
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*'
_output_shapes
:���������*$
_read_only_resource_inputs
**
config_proto

GPU 

CPU2J 8*O
fJRH
F__inference_dense_880_layer_call_and_return_conditional_losses_10167222
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������(::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������(
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: 
�'
�
J__inference_sequential_87_layer_call_and_return_conditional_losses_1016892

inputs
dense_874_1016856
dense_874_1016858
dense_875_1016861
dense_875_1016863
dense_876_1016866
dense_876_1016868
dense_877_1016871
dense_877_1016873
dense_878_1016876
dense_878_1016878
dense_879_1016881
dense_879_1016883
dense_880_1016886
dense_880_1016888
identity��!dense_874/StatefulPartitionedCall�!dense_875/StatefulPartitionedCall�!dense_876/StatefulPartitionedCall�!dense_877/StatefulPartitionedCall�!dense_878/StatefulPartitionedCall�!dense_879/StatefulPartitionedCall�!dense_880/StatefulPartitionedCall�
!dense_874/StatefulPartitionedCallStatefulPartitionedCallinputsdense_874_1016856dense_874_1016858*
Tin
2*
Tout
2*'
_output_shapes
:���������(*$
_read_only_resource_inputs
**
config_proto

GPU 

CPU2J 8*O
fJRH
F__inference_dense_874_layer_call_and_return_conditional_losses_10165612#
!dense_874/StatefulPartitionedCall�
!dense_875/StatefulPartitionedCallStatefulPartitionedCall*dense_874/StatefulPartitionedCall:output:0dense_875_1016861dense_875_1016863*
Tin
2*
Tout
2*'
_output_shapes
:���������(*$
_read_only_resource_inputs
**
config_proto

GPU 

CPU2J 8*O
fJRH
F__inference_dense_875_layer_call_and_return_conditional_losses_10165882#
!dense_875/StatefulPartitionedCall�
!dense_876/StatefulPartitionedCallStatefulPartitionedCall*dense_875/StatefulPartitionedCall:output:0dense_876_1016866dense_876_1016868*
Tin
2*
Tout
2*'
_output_shapes
:���������(*$
_read_only_resource_inputs
**
config_proto

GPU 

CPU2J 8*O
fJRH
F__inference_dense_876_layer_call_and_return_conditional_losses_10166152#
!dense_876/StatefulPartitionedCall�
!dense_877/StatefulPartitionedCallStatefulPartitionedCall*dense_876/StatefulPartitionedCall:output:0dense_877_1016871dense_877_1016873*
Tin
2*
Tout
2*'
_output_shapes
:���������(*$
_read_only_resource_inputs
**
config_proto

GPU 

CPU2J 8*O
fJRH
F__inference_dense_877_layer_call_and_return_conditional_losses_10166422#
!dense_877/StatefulPartitionedCall�
!dense_878/StatefulPartitionedCallStatefulPartitionedCall*dense_877/StatefulPartitionedCall:output:0dense_878_1016876dense_878_1016878*
Tin
2*
Tout
2*'
_output_shapes
:���������(*$
_read_only_resource_inputs
**
config_proto

GPU 

CPU2J 8*O
fJRH
F__inference_dense_878_layer_call_and_return_conditional_losses_10166692#
!dense_878/StatefulPartitionedCall�
!dense_879/StatefulPartitionedCallStatefulPartitionedCall*dense_878/StatefulPartitionedCall:output:0dense_879_1016881dense_879_1016883*
Tin
2*
Tout
2*'
_output_shapes
:���������(*$
_read_only_resource_inputs
**
config_proto

GPU 

CPU2J 8*O
fJRH
F__inference_dense_879_layer_call_and_return_conditional_losses_10166962#
!dense_879/StatefulPartitionedCall�
!dense_880/StatefulPartitionedCallStatefulPartitionedCall*dense_879/StatefulPartitionedCall:output:0dense_880_1016886dense_880_1016888*
Tin
2*
Tout
2*'
_output_shapes
:���������*$
_read_only_resource_inputs
**
config_proto

GPU 

CPU2J 8*O
fJRH
F__inference_dense_880_layer_call_and_return_conditional_losses_10167222#
!dense_880/StatefulPartitionedCall�
IdentityIdentity*dense_880/StatefulPartitionedCall:output:0"^dense_874/StatefulPartitionedCall"^dense_875/StatefulPartitionedCall"^dense_876/StatefulPartitionedCall"^dense_877/StatefulPartitionedCall"^dense_878/StatefulPartitionedCall"^dense_879/StatefulPartitionedCall"^dense_880/StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*^
_input_shapesM
K:���������::::::::::::::2F
!dense_874/StatefulPartitionedCall!dense_874/StatefulPartitionedCall2F
!dense_875/StatefulPartitionedCall!dense_875/StatefulPartitionedCall2F
!dense_876/StatefulPartitionedCall!dense_876/StatefulPartitionedCall2F
!dense_877/StatefulPartitionedCall!dense_877/StatefulPartitionedCall2F
!dense_878/StatefulPartitionedCall!dense_878/StatefulPartitionedCall2F
!dense_879/StatefulPartitionedCall!dense_879/StatefulPartitionedCall2F
!dense_880/StatefulPartitionedCall!dense_880/StatefulPartitionedCall:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :	

_output_shapes
: :


_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: 
�t
�
 __inference__traced_save_1017447
file_prefix/
+savev2_dense_874_kernel_read_readvariableop-
)savev2_dense_874_bias_read_readvariableop/
+savev2_dense_875_kernel_read_readvariableop-
)savev2_dense_875_bias_read_readvariableop/
+savev2_dense_876_kernel_read_readvariableop-
)savev2_dense_876_bias_read_readvariableop/
+savev2_dense_877_kernel_read_readvariableop-
)savev2_dense_877_bias_read_readvariableop/
+savev2_dense_878_kernel_read_readvariableop-
)savev2_dense_878_bias_read_readvariableop/
+savev2_dense_879_kernel_read_readvariableop-
)savev2_dense_879_bias_read_readvariableop/
+savev2_dense_880_kernel_read_readvariableop-
)savev2_dense_880_bias_read_readvariableop,
(savev2_adadelta_iter_read_readvariableop	-
)savev2_adadelta_decay_read_readvariableop5
1savev2_adadelta_learning_rate_read_readvariableop+
'savev2_adadelta_rho_read_readvariableop$
 savev2_total_read_readvariableop$
 savev2_count_read_readvariableop&
"savev2_total_1_read_readvariableop&
"savev2_count_1_read_readvariableopC
?savev2_adadelta_dense_874_kernel_accum_grad_read_readvariableopA
=savev2_adadelta_dense_874_bias_accum_grad_read_readvariableopC
?savev2_adadelta_dense_875_kernel_accum_grad_read_readvariableopA
=savev2_adadelta_dense_875_bias_accum_grad_read_readvariableopC
?savev2_adadelta_dense_876_kernel_accum_grad_read_readvariableopA
=savev2_adadelta_dense_876_bias_accum_grad_read_readvariableopC
?savev2_adadelta_dense_877_kernel_accum_grad_read_readvariableopA
=savev2_adadelta_dense_877_bias_accum_grad_read_readvariableopC
?savev2_adadelta_dense_878_kernel_accum_grad_read_readvariableopA
=savev2_adadelta_dense_878_bias_accum_grad_read_readvariableopC
?savev2_adadelta_dense_879_kernel_accum_grad_read_readvariableopA
=savev2_adadelta_dense_879_bias_accum_grad_read_readvariableopC
?savev2_adadelta_dense_880_kernel_accum_grad_read_readvariableopA
=savev2_adadelta_dense_880_bias_accum_grad_read_readvariableopB
>savev2_adadelta_dense_874_kernel_accum_var_read_readvariableop@
<savev2_adadelta_dense_874_bias_accum_var_read_readvariableopB
>savev2_adadelta_dense_875_kernel_accum_var_read_readvariableop@
<savev2_adadelta_dense_875_bias_accum_var_read_readvariableopB
>savev2_adadelta_dense_876_kernel_accum_var_read_readvariableop@
<savev2_adadelta_dense_876_bias_accum_var_read_readvariableopB
>savev2_adadelta_dense_877_kernel_accum_var_read_readvariableop@
<savev2_adadelta_dense_877_bias_accum_var_read_readvariableopB
>savev2_adadelta_dense_878_kernel_accum_var_read_readvariableop@
<savev2_adadelta_dense_878_bias_accum_var_read_readvariableopB
>savev2_adadelta_dense_879_kernel_accum_var_read_readvariableop@
<savev2_adadelta_dense_879_bias_accum_var_read_readvariableopB
>savev2_adadelta_dense_880_kernel_accum_var_read_readvariableop@
<savev2_adadelta_dense_880_bias_accum_var_read_readvariableop
savev2_1_const

identity_1��MergeV2Checkpoints�SaveV2�SaveV2_1�
StaticRegexFullMatchStaticRegexFullMatchfile_prefix"/device:CPU:**
_output_shapes
: *
pattern
^s3://.*2
StaticRegexFullMatchc
ConstConst"/device:CPU:**
_output_shapes
: *
dtype0*
valueB B.part2
Const�
Const_1Const"/device:CPU:**
_output_shapes
: *
dtype0*<
value3B1 B+_temp_f5d57b28199e4a75a7c738eaa85fc243/part2	
Const_1�
SelectSelectStaticRegexFullMatch:output:0Const:output:0Const_1:output:0"/device:CPU:**
T0*
_output_shapes
: 2
Selectt

StringJoin
StringJoinfile_prefixSelect:output:0"/device:CPU:**
N*
_output_shapes
: 2

StringJoinZ

num_shardsConst*
_output_shapes
: *
dtype0*
value	B :2

num_shards
ShardedFilename/shardConst"/device:CPU:0*
_output_shapes
: *
dtype0*
value	B : 2
ShardedFilename/shard�
ShardedFilenameShardedFilenameStringJoin:output:0ShardedFilename/shard:output:0num_shards:output:0"/device:CPU:0*
_output_shapes
: 2
ShardedFilename�
SaveV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:2*
dtype0*�
value�B�2B6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-2/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-3/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-3/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-4/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-4/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-5/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-5/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-6/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-6/bias/.ATTRIBUTES/VARIABLE_VALUEB)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUEB*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUEB(optimizer/rho/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUEB[layer_with_weights-0/kernel/.OPTIMIZER_SLOT/optimizer/accum_grad/.ATTRIBUTES/VARIABLE_VALUEBYlayer_with_weights-0/bias/.OPTIMIZER_SLOT/optimizer/accum_grad/.ATTRIBUTES/VARIABLE_VALUEB[layer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/accum_grad/.ATTRIBUTES/VARIABLE_VALUEBYlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/accum_grad/.ATTRIBUTES/VARIABLE_VALUEB[layer_with_weights-2/kernel/.OPTIMIZER_SLOT/optimizer/accum_grad/.ATTRIBUTES/VARIABLE_VALUEBYlayer_with_weights-2/bias/.OPTIMIZER_SLOT/optimizer/accum_grad/.ATTRIBUTES/VARIABLE_VALUEB[layer_with_weights-3/kernel/.OPTIMIZER_SLOT/optimizer/accum_grad/.ATTRIBUTES/VARIABLE_VALUEBYlayer_with_weights-3/bias/.OPTIMIZER_SLOT/optimizer/accum_grad/.ATTRIBUTES/VARIABLE_VALUEB[layer_with_weights-4/kernel/.OPTIMIZER_SLOT/optimizer/accum_grad/.ATTRIBUTES/VARIABLE_VALUEBYlayer_with_weights-4/bias/.OPTIMIZER_SLOT/optimizer/accum_grad/.ATTRIBUTES/VARIABLE_VALUEB[layer_with_weights-5/kernel/.OPTIMIZER_SLOT/optimizer/accum_grad/.ATTRIBUTES/VARIABLE_VALUEBYlayer_with_weights-5/bias/.OPTIMIZER_SLOT/optimizer/accum_grad/.ATTRIBUTES/VARIABLE_VALUEB[layer_with_weights-6/kernel/.OPTIMIZER_SLOT/optimizer/accum_grad/.ATTRIBUTES/VARIABLE_VALUEBYlayer_with_weights-6/bias/.OPTIMIZER_SLOT/optimizer/accum_grad/.ATTRIBUTES/VARIABLE_VALUEBZlayer_with_weights-0/kernel/.OPTIMIZER_SLOT/optimizer/accum_var/.ATTRIBUTES/VARIABLE_VALUEBXlayer_with_weights-0/bias/.OPTIMIZER_SLOT/optimizer/accum_var/.ATTRIBUTES/VARIABLE_VALUEBZlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/accum_var/.ATTRIBUTES/VARIABLE_VALUEBXlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/accum_var/.ATTRIBUTES/VARIABLE_VALUEBZlayer_with_weights-2/kernel/.OPTIMIZER_SLOT/optimizer/accum_var/.ATTRIBUTES/VARIABLE_VALUEBXlayer_with_weights-2/bias/.OPTIMIZER_SLOT/optimizer/accum_var/.ATTRIBUTES/VARIABLE_VALUEBZlayer_with_weights-3/kernel/.OPTIMIZER_SLOT/optimizer/accum_var/.ATTRIBUTES/VARIABLE_VALUEBXlayer_with_weights-3/bias/.OPTIMIZER_SLOT/optimizer/accum_var/.ATTRIBUTES/VARIABLE_VALUEBZlayer_with_weights-4/kernel/.OPTIMIZER_SLOT/optimizer/accum_var/.ATTRIBUTES/VARIABLE_VALUEBXlayer_with_weights-4/bias/.OPTIMIZER_SLOT/optimizer/accum_var/.ATTRIBUTES/VARIABLE_VALUEBZlayer_with_weights-5/kernel/.OPTIMIZER_SLOT/optimizer/accum_var/.ATTRIBUTES/VARIABLE_VALUEBXlayer_with_weights-5/bias/.OPTIMIZER_SLOT/optimizer/accum_var/.ATTRIBUTES/VARIABLE_VALUEBZlayer_with_weights-6/kernel/.OPTIMIZER_SLOT/optimizer/accum_var/.ATTRIBUTES/VARIABLE_VALUEBXlayer_with_weights-6/bias/.OPTIMIZER_SLOT/optimizer/accum_var/.ATTRIBUTES/VARIABLE_VALUE2
SaveV2/tensor_names�
SaveV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:2*
dtype0*w
valuenBl2B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B 2
SaveV2/shape_and_slices�
SaveV2SaveV2ShardedFilename:filename:0SaveV2/tensor_names:output:0 SaveV2/shape_and_slices:output:0+savev2_dense_874_kernel_read_readvariableop)savev2_dense_874_bias_read_readvariableop+savev2_dense_875_kernel_read_readvariableop)savev2_dense_875_bias_read_readvariableop+savev2_dense_876_kernel_read_readvariableop)savev2_dense_876_bias_read_readvariableop+savev2_dense_877_kernel_read_readvariableop)savev2_dense_877_bias_read_readvariableop+savev2_dense_878_kernel_read_readvariableop)savev2_dense_878_bias_read_readvariableop+savev2_dense_879_kernel_read_readvariableop)savev2_dense_879_bias_read_readvariableop+savev2_dense_880_kernel_read_readvariableop)savev2_dense_880_bias_read_readvariableop(savev2_adadelta_iter_read_readvariableop)savev2_adadelta_decay_read_readvariableop1savev2_adadelta_learning_rate_read_readvariableop'savev2_adadelta_rho_read_readvariableop savev2_total_read_readvariableop savev2_count_read_readvariableop"savev2_total_1_read_readvariableop"savev2_count_1_read_readvariableop?savev2_adadelta_dense_874_kernel_accum_grad_read_readvariableop=savev2_adadelta_dense_874_bias_accum_grad_read_readvariableop?savev2_adadelta_dense_875_kernel_accum_grad_read_readvariableop=savev2_adadelta_dense_875_bias_accum_grad_read_readvariableop?savev2_adadelta_dense_876_kernel_accum_grad_read_readvariableop=savev2_adadelta_dense_876_bias_accum_grad_read_readvariableop?savev2_adadelta_dense_877_kernel_accum_grad_read_readvariableop=savev2_adadelta_dense_877_bias_accum_grad_read_readvariableop?savev2_adadelta_dense_878_kernel_accum_grad_read_readvariableop=savev2_adadelta_dense_878_bias_accum_grad_read_readvariableop?savev2_adadelta_dense_879_kernel_accum_grad_read_readvariableop=savev2_adadelta_dense_879_bias_accum_grad_read_readvariableop?savev2_adadelta_dense_880_kernel_accum_grad_read_readvariableop=savev2_adadelta_dense_880_bias_accum_grad_read_readvariableop>savev2_adadelta_dense_874_kernel_accum_var_read_readvariableop<savev2_adadelta_dense_874_bias_accum_var_read_readvariableop>savev2_adadelta_dense_875_kernel_accum_var_read_readvariableop<savev2_adadelta_dense_875_bias_accum_var_read_readvariableop>savev2_adadelta_dense_876_kernel_accum_var_read_readvariableop<savev2_adadelta_dense_876_bias_accum_var_read_readvariableop>savev2_adadelta_dense_877_kernel_accum_var_read_readvariableop<savev2_adadelta_dense_877_bias_accum_var_read_readvariableop>savev2_adadelta_dense_878_kernel_accum_var_read_readvariableop<savev2_adadelta_dense_878_bias_accum_var_read_readvariableop>savev2_adadelta_dense_879_kernel_accum_var_read_readvariableop<savev2_adadelta_dense_879_bias_accum_var_read_readvariableop>savev2_adadelta_dense_880_kernel_accum_var_read_readvariableop<savev2_adadelta_dense_880_bias_accum_var_read_readvariableop"/device:CPU:0*
_output_shapes
 *@
dtypes6
422	2
SaveV2�
ShardedFilename_1/shardConst"/device:CPU:0*
_output_shapes
: *
dtype0*
value	B :2
ShardedFilename_1/shard�
ShardedFilename_1ShardedFilenameStringJoin:output:0 ShardedFilename_1/shard:output:0num_shards:output:0"/device:CPU:0*
_output_shapes
: 2
ShardedFilename_1�
SaveV2_1/tensor_namesConst"/device:CPU:0*
_output_shapes
:*
dtype0*1
value(B&B_CHECKPOINTABLE_OBJECT_GRAPH2
SaveV2_1/tensor_names�
SaveV2_1/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:*
dtype0*
valueB
B 2
SaveV2_1/shape_and_slices�
SaveV2_1SaveV2ShardedFilename_1:filename:0SaveV2_1/tensor_names:output:0"SaveV2_1/shape_and_slices:output:0savev2_1_const^SaveV2"/device:CPU:0*
_output_shapes
 *
dtypes
22

SaveV2_1�
&MergeV2Checkpoints/checkpoint_prefixesPackShardedFilename:filename:0ShardedFilename_1:filename:0^SaveV2	^SaveV2_1"/device:CPU:0*
N*
T0*
_output_shapes
:2(
&MergeV2Checkpoints/checkpoint_prefixes�
MergeV2CheckpointsMergeV2Checkpoints/MergeV2Checkpoints/checkpoint_prefixes:output:0file_prefix	^SaveV2_1"/device:CPU:0*
_output_shapes
 2
MergeV2Checkpointsr
IdentityIdentityfile_prefix^MergeV2Checkpoints"/device:CPU:0*
T0*
_output_shapes
: 2

Identity�

Identity_1IdentityIdentity:output:0^MergeV2Checkpoints^SaveV2	^SaveV2_1*
T0*
_output_shapes
: 2

Identity_1"!

identity_1Identity_1:output:0*�
_input_shapes�
�: :(:(:((:(:((:(:((:(:((:(:((:(:(:: : : : : : : : :(:(:((:(:((:(:((:(:((:(:((:(:(::(:(:((:(:((:(:((:(:((:(:((:(:(:: 2(
MergeV2CheckpointsMergeV2Checkpoints2
SaveV2SaveV22
SaveV2_1SaveV2_1:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix:$ 

_output_shapes

:(: 

_output_shapes
:(:$ 

_output_shapes

:((: 

_output_shapes
:(:$ 

_output_shapes

:((: 

_output_shapes
:(:$ 

_output_shapes

:((: 

_output_shapes
:(:$	 

_output_shapes

:((: 


_output_shapes
:(:$ 

_output_shapes

:((: 

_output_shapes
:(:$ 

_output_shapes

:(: 

_output_shapes
::

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :$ 

_output_shapes

:(: 

_output_shapes
:(:$ 

_output_shapes

:((: 

_output_shapes
:(:$ 

_output_shapes

:((: 

_output_shapes
:(:$ 

_output_shapes

:((: 

_output_shapes
:(:$ 

_output_shapes

:((:  

_output_shapes
:(:$! 

_output_shapes

:((: "

_output_shapes
:(:$# 

_output_shapes

:(: $

_output_shapes
::$% 

_output_shapes

:(: &

_output_shapes
:(:$' 

_output_shapes

:((: (

_output_shapes
:(:$) 

_output_shapes

:((: *

_output_shapes
:(:$+ 

_output_shapes

:((: ,

_output_shapes
:(:$- 

_output_shapes

:((: .

_output_shapes
:(:$/ 

_output_shapes

:((: 0

_output_shapes
:(:$1 

_output_shapes

:(: 2

_output_shapes
::3

_output_shapes
: 
�
�
F__inference_dense_878_layer_call_and_return_conditional_losses_1016669

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity��
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:((*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2
MatMul�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:(*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������(2
Reluf
IdentityIdentityRelu:activations:0*
T0*'
_output_shapes
:���������(2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������(:::O K
'
_output_shapes
:���������(
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: 
�
�
+__inference_dense_876_layer_call_fn_1017191

inputs
unknown
	unknown_0
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*'
_output_shapes
:���������(*$
_read_only_resource_inputs
**
config_proto

GPU 

CPU2J 8*O
fJRH
F__inference_dense_876_layer_call_and_return_conditional_losses_10166152
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������(2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������(::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������(
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: 
�
�
F__inference_dense_875_layer_call_and_return_conditional_losses_1016588

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity��
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:((*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2
MatMul�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:(*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������(2
Reluf
IdentityIdentityRelu:activations:0*
T0*'
_output_shapes
:���������(2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������(:::O K
'
_output_shapes
:���������(
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: 
�'
�
J__inference_sequential_87_layer_call_and_return_conditional_losses_1016820

inputs
dense_874_1016784
dense_874_1016786
dense_875_1016789
dense_875_1016791
dense_876_1016794
dense_876_1016796
dense_877_1016799
dense_877_1016801
dense_878_1016804
dense_878_1016806
dense_879_1016809
dense_879_1016811
dense_880_1016814
dense_880_1016816
identity��!dense_874/StatefulPartitionedCall�!dense_875/StatefulPartitionedCall�!dense_876/StatefulPartitionedCall�!dense_877/StatefulPartitionedCall�!dense_878/StatefulPartitionedCall�!dense_879/StatefulPartitionedCall�!dense_880/StatefulPartitionedCall�
!dense_874/StatefulPartitionedCallStatefulPartitionedCallinputsdense_874_1016784dense_874_1016786*
Tin
2*
Tout
2*'
_output_shapes
:���������(*$
_read_only_resource_inputs
**
config_proto

GPU 

CPU2J 8*O
fJRH
F__inference_dense_874_layer_call_and_return_conditional_losses_10165612#
!dense_874/StatefulPartitionedCall�
!dense_875/StatefulPartitionedCallStatefulPartitionedCall*dense_874/StatefulPartitionedCall:output:0dense_875_1016789dense_875_1016791*
Tin
2*
Tout
2*'
_output_shapes
:���������(*$
_read_only_resource_inputs
**
config_proto

GPU 

CPU2J 8*O
fJRH
F__inference_dense_875_layer_call_and_return_conditional_losses_10165882#
!dense_875/StatefulPartitionedCall�
!dense_876/StatefulPartitionedCallStatefulPartitionedCall*dense_875/StatefulPartitionedCall:output:0dense_876_1016794dense_876_1016796*
Tin
2*
Tout
2*'
_output_shapes
:���������(*$
_read_only_resource_inputs
**
config_proto

GPU 

CPU2J 8*O
fJRH
F__inference_dense_876_layer_call_and_return_conditional_losses_10166152#
!dense_876/StatefulPartitionedCall�
!dense_877/StatefulPartitionedCallStatefulPartitionedCall*dense_876/StatefulPartitionedCall:output:0dense_877_1016799dense_877_1016801*
Tin
2*
Tout
2*'
_output_shapes
:���������(*$
_read_only_resource_inputs
**
config_proto

GPU 

CPU2J 8*O
fJRH
F__inference_dense_877_layer_call_and_return_conditional_losses_10166422#
!dense_877/StatefulPartitionedCall�
!dense_878/StatefulPartitionedCallStatefulPartitionedCall*dense_877/StatefulPartitionedCall:output:0dense_878_1016804dense_878_1016806*
Tin
2*
Tout
2*'
_output_shapes
:���������(*$
_read_only_resource_inputs
**
config_proto

GPU 

CPU2J 8*O
fJRH
F__inference_dense_878_layer_call_and_return_conditional_losses_10166692#
!dense_878/StatefulPartitionedCall�
!dense_879/StatefulPartitionedCallStatefulPartitionedCall*dense_878/StatefulPartitionedCall:output:0dense_879_1016809dense_879_1016811*
Tin
2*
Tout
2*'
_output_shapes
:���������(*$
_read_only_resource_inputs
**
config_proto

GPU 

CPU2J 8*O
fJRH
F__inference_dense_879_layer_call_and_return_conditional_losses_10166962#
!dense_879/StatefulPartitionedCall�
!dense_880/StatefulPartitionedCallStatefulPartitionedCall*dense_879/StatefulPartitionedCall:output:0dense_880_1016814dense_880_1016816*
Tin
2*
Tout
2*'
_output_shapes
:���������*$
_read_only_resource_inputs
**
config_proto

GPU 

CPU2J 8*O
fJRH
F__inference_dense_880_layer_call_and_return_conditional_losses_10167222#
!dense_880/StatefulPartitionedCall�
IdentityIdentity*dense_880/StatefulPartitionedCall:output:0"^dense_874/StatefulPartitionedCall"^dense_875/StatefulPartitionedCall"^dense_876/StatefulPartitionedCall"^dense_877/StatefulPartitionedCall"^dense_878/StatefulPartitionedCall"^dense_879/StatefulPartitionedCall"^dense_880/StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*^
_input_shapesM
K:���������::::::::::::::2F
!dense_874/StatefulPartitionedCall!dense_874/StatefulPartitionedCall2F
!dense_875/StatefulPartitionedCall!dense_875/StatefulPartitionedCall2F
!dense_876/StatefulPartitionedCall!dense_876/StatefulPartitionedCall2F
!dense_877/StatefulPartitionedCall!dense_877/StatefulPartitionedCall2F
!dense_878/StatefulPartitionedCall!dense_878/StatefulPartitionedCall2F
!dense_879/StatefulPartitionedCall!dense_879/StatefulPartitionedCall2F
!dense_880/StatefulPartitionedCall!dense_880/StatefulPartitionedCall:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :	

_output_shapes
: :


_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: 
�
�
/__inference_sequential_87_layer_call_fn_1016851
dense_874_input
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
	unknown_9

unknown_10

unknown_11

unknown_12
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCalldense_874_inputunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12*
Tin
2*
Tout
2*'
_output_shapes
:���������*0
_read_only_resource_inputs
	
**
config_proto

GPU 

CPU2J 8*S
fNRL
J__inference_sequential_87_layer_call_and_return_conditional_losses_10168202
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*^
_input_shapesM
K:���������::::::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:X T
'
_output_shapes
:���������
)
_user_specified_namedense_874_input:

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :	

_output_shapes
: :


_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: 
�
�
F__inference_dense_876_layer_call_and_return_conditional_losses_1017182

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity��
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:((*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2
MatMul�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:(*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������(2
Reluf
IdentityIdentityRelu:activations:0*
T0*'
_output_shapes
:���������(2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������(:::O K
'
_output_shapes
:���������(
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: 
�
�
F__inference_dense_876_layer_call_and_return_conditional_losses_1016615

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity��
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:((*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2
MatMul�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:(*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������(2
Reluf
IdentityIdentityRelu:activations:0*
T0*'
_output_shapes
:���������(2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������(:::O K
'
_output_shapes
:���������(
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: 
�
�
/__inference_sequential_87_layer_call_fn_1017132

inputs
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
	unknown_9

unknown_10

unknown_11

unknown_12
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12*
Tin
2*
Tout
2*'
_output_shapes
:���������*0
_read_only_resource_inputs
	
**
config_proto

GPU 

CPU2J 8*S
fNRL
J__inference_sequential_87_layer_call_and_return_conditional_losses_10168922
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*^
_input_shapesM
K:���������::::::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :	

_output_shapes
: :


_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: 
�C
�
"__inference__wrapped_model_1016547
dense_874_input:
6sequential_87_dense_874_matmul_readvariableop_resource;
7sequential_87_dense_874_biasadd_readvariableop_resource:
6sequential_87_dense_875_matmul_readvariableop_resource;
7sequential_87_dense_875_biasadd_readvariableop_resource:
6sequential_87_dense_876_matmul_readvariableop_resource;
7sequential_87_dense_876_biasadd_readvariableop_resource:
6sequential_87_dense_877_matmul_readvariableop_resource;
7sequential_87_dense_877_biasadd_readvariableop_resource:
6sequential_87_dense_878_matmul_readvariableop_resource;
7sequential_87_dense_878_biasadd_readvariableop_resource:
6sequential_87_dense_879_matmul_readvariableop_resource;
7sequential_87_dense_879_biasadd_readvariableop_resource:
6sequential_87_dense_880_matmul_readvariableop_resource;
7sequential_87_dense_880_biasadd_readvariableop_resource
identity��
-sequential_87/dense_874/MatMul/ReadVariableOpReadVariableOp6sequential_87_dense_874_matmul_readvariableop_resource*
_output_shapes

:(*
dtype02/
-sequential_87/dense_874/MatMul/ReadVariableOp�
sequential_87/dense_874/MatMulMatMuldense_874_input5sequential_87/dense_874/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2 
sequential_87/dense_874/MatMul�
.sequential_87/dense_874/BiasAdd/ReadVariableOpReadVariableOp7sequential_87_dense_874_biasadd_readvariableop_resource*
_output_shapes
:(*
dtype020
.sequential_87/dense_874/BiasAdd/ReadVariableOp�
sequential_87/dense_874/BiasAddBiasAdd(sequential_87/dense_874/MatMul:product:06sequential_87/dense_874/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2!
sequential_87/dense_874/BiasAdd�
-sequential_87/dense_875/MatMul/ReadVariableOpReadVariableOp6sequential_87_dense_875_matmul_readvariableop_resource*
_output_shapes

:((*
dtype02/
-sequential_87/dense_875/MatMul/ReadVariableOp�
sequential_87/dense_875/MatMulMatMul(sequential_87/dense_874/BiasAdd:output:05sequential_87/dense_875/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2 
sequential_87/dense_875/MatMul�
.sequential_87/dense_875/BiasAdd/ReadVariableOpReadVariableOp7sequential_87_dense_875_biasadd_readvariableop_resource*
_output_shapes
:(*
dtype020
.sequential_87/dense_875/BiasAdd/ReadVariableOp�
sequential_87/dense_875/BiasAddBiasAdd(sequential_87/dense_875/MatMul:product:06sequential_87/dense_875/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2!
sequential_87/dense_875/BiasAdd�
sequential_87/dense_875/ReluRelu(sequential_87/dense_875/BiasAdd:output:0*
T0*'
_output_shapes
:���������(2
sequential_87/dense_875/Relu�
-sequential_87/dense_876/MatMul/ReadVariableOpReadVariableOp6sequential_87_dense_876_matmul_readvariableop_resource*
_output_shapes

:((*
dtype02/
-sequential_87/dense_876/MatMul/ReadVariableOp�
sequential_87/dense_876/MatMulMatMul*sequential_87/dense_875/Relu:activations:05sequential_87/dense_876/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2 
sequential_87/dense_876/MatMul�
.sequential_87/dense_876/BiasAdd/ReadVariableOpReadVariableOp7sequential_87_dense_876_biasadd_readvariableop_resource*
_output_shapes
:(*
dtype020
.sequential_87/dense_876/BiasAdd/ReadVariableOp�
sequential_87/dense_876/BiasAddBiasAdd(sequential_87/dense_876/MatMul:product:06sequential_87/dense_876/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2!
sequential_87/dense_876/BiasAdd�
sequential_87/dense_876/ReluRelu(sequential_87/dense_876/BiasAdd:output:0*
T0*'
_output_shapes
:���������(2
sequential_87/dense_876/Relu�
-sequential_87/dense_877/MatMul/ReadVariableOpReadVariableOp6sequential_87_dense_877_matmul_readvariableop_resource*
_output_shapes

:((*
dtype02/
-sequential_87/dense_877/MatMul/ReadVariableOp�
sequential_87/dense_877/MatMulMatMul*sequential_87/dense_876/Relu:activations:05sequential_87/dense_877/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2 
sequential_87/dense_877/MatMul�
.sequential_87/dense_877/BiasAdd/ReadVariableOpReadVariableOp7sequential_87_dense_877_biasadd_readvariableop_resource*
_output_shapes
:(*
dtype020
.sequential_87/dense_877/BiasAdd/ReadVariableOp�
sequential_87/dense_877/BiasAddBiasAdd(sequential_87/dense_877/MatMul:product:06sequential_87/dense_877/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2!
sequential_87/dense_877/BiasAdd�
sequential_87/dense_877/ReluRelu(sequential_87/dense_877/BiasAdd:output:0*
T0*'
_output_shapes
:���������(2
sequential_87/dense_877/Relu�
-sequential_87/dense_878/MatMul/ReadVariableOpReadVariableOp6sequential_87_dense_878_matmul_readvariableop_resource*
_output_shapes

:((*
dtype02/
-sequential_87/dense_878/MatMul/ReadVariableOp�
sequential_87/dense_878/MatMulMatMul*sequential_87/dense_877/Relu:activations:05sequential_87/dense_878/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2 
sequential_87/dense_878/MatMul�
.sequential_87/dense_878/BiasAdd/ReadVariableOpReadVariableOp7sequential_87_dense_878_biasadd_readvariableop_resource*
_output_shapes
:(*
dtype020
.sequential_87/dense_878/BiasAdd/ReadVariableOp�
sequential_87/dense_878/BiasAddBiasAdd(sequential_87/dense_878/MatMul:product:06sequential_87/dense_878/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2!
sequential_87/dense_878/BiasAdd�
sequential_87/dense_878/ReluRelu(sequential_87/dense_878/BiasAdd:output:0*
T0*'
_output_shapes
:���������(2
sequential_87/dense_878/Relu�
-sequential_87/dense_879/MatMul/ReadVariableOpReadVariableOp6sequential_87_dense_879_matmul_readvariableop_resource*
_output_shapes

:((*
dtype02/
-sequential_87/dense_879/MatMul/ReadVariableOp�
sequential_87/dense_879/MatMulMatMul*sequential_87/dense_878/Relu:activations:05sequential_87/dense_879/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2 
sequential_87/dense_879/MatMul�
.sequential_87/dense_879/BiasAdd/ReadVariableOpReadVariableOp7sequential_87_dense_879_biasadd_readvariableop_resource*
_output_shapes
:(*
dtype020
.sequential_87/dense_879/BiasAdd/ReadVariableOp�
sequential_87/dense_879/BiasAddBiasAdd(sequential_87/dense_879/MatMul:product:06sequential_87/dense_879/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������(2!
sequential_87/dense_879/BiasAdd�
sequential_87/dense_879/ReluRelu(sequential_87/dense_879/BiasAdd:output:0*
T0*'
_output_shapes
:���������(2
sequential_87/dense_879/Relu�
-sequential_87/dense_880/MatMul/ReadVariableOpReadVariableOp6sequential_87_dense_880_matmul_readvariableop_resource*
_output_shapes

:(*
dtype02/
-sequential_87/dense_880/MatMul/ReadVariableOp�
sequential_87/dense_880/MatMulMatMul*sequential_87/dense_879/Relu:activations:05sequential_87/dense_880/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2 
sequential_87/dense_880/MatMul�
.sequential_87/dense_880/BiasAdd/ReadVariableOpReadVariableOp7sequential_87_dense_880_biasadd_readvariableop_resource*
_output_shapes
:*
dtype020
.sequential_87/dense_880/BiasAdd/ReadVariableOp�
sequential_87/dense_880/BiasAddBiasAdd(sequential_87/dense_880/MatMul:product:06sequential_87/dense_880/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2!
sequential_87/dense_880/BiasAdd|
IdentityIdentity(sequential_87/dense_880/BiasAdd:output:0*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*^
_input_shapesM
K:���������:::::::::::::::X T
'
_output_shapes
:���������
)
_user_specified_namedense_874_input:

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :	

_output_shapes
: :


_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: "�L
saver_filename:0StatefulPartitionedCall_1:0StatefulPartitionedCall_28"
saved_model_main_op

NoOp*>
__saved_model_init_op%#
__saved_model_init_op

NoOp*�
serving_default�
K
dense_874_input8
!serving_default_dense_874_input:0���������=
	dense_8800
StatefulPartitionedCall:0���������tensorflow/serving/predict:��
�;
layer_with_weights-0
layer-0
layer_with_weights-1
layer-1
layer_with_weights-2
layer-2
layer_with_weights-3
layer-3
layer_with_weights-4
layer-4
layer_with_weights-5
layer-5
layer_with_weights-6
layer-6
	optimizer
	regularization_losses

	variables
trainable_variables
	keras_api

signatures
�__call__
+�&call_and_return_all_conditional_losses
�_default_save_signature"�8
_tf_keras_sequential�7{"class_name": "Sequential", "name": "sequential_87", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "config": {"name": "sequential_87", "layers": [{"class_name": "Dense", "config": {"name": "dense_874", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 1]}, "dtype": "float32", "units": 40, "activation": "linear", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}, {"class_name": "Dense", "config": {"name": "dense_875", "trainable": true, "dtype": "float32", "units": 40, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}, {"class_name": "Dense", "config": {"name": "dense_876", "trainable": true, "dtype": "float32", "units": 40, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}, {"class_name": "Dense", "config": {"name": "dense_877", "trainable": true, "dtype": "float32", "units": 40, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}, {"class_name": "Dense", "config": {"name": "dense_878", "trainable": true, "dtype": "float32", "units": 40, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}, {"class_name": "Dense", "config": {"name": "dense_879", "trainable": true, "dtype": "float32", "units": 40, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}, {"class_name": "Dense", "config": {"name": "dense_880", "trainable": true, "dtype": "float32", "units": 1, "activation": "linear", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}], "build_input_shape": {"class_name": "TensorShape", "items": [null, 1]}}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 1}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 1]}, "is_graph_network": true, "keras_version": "2.3.0-tf", "backend": "tensorflow", "model_config": {"class_name": "Sequential", "config": {"name": "sequential_87", "layers": [{"class_name": "Dense", "config": {"name": "dense_874", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 1]}, "dtype": "float32", "units": 40, "activation": "linear", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}, {"class_name": "Dense", "config": {"name": "dense_875", "trainable": true, "dtype": "float32", "units": 40, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}, {"class_name": "Dense", "config": {"name": "dense_876", "trainable": true, "dtype": "float32", "units": 40, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}, {"class_name": "Dense", "config": {"name": "dense_877", "trainable": true, "dtype": "float32", "units": 40, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}, {"class_name": "Dense", "config": {"name": "dense_878", "trainable": true, "dtype": "float32", "units": 40, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}, {"class_name": "Dense", "config": {"name": "dense_879", "trainable": true, "dtype": "float32", "units": 40, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}, {"class_name": "Dense", "config": {"name": "dense_880", "trainable": true, "dtype": "float32", "units": 1, "activation": "linear", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}], "build_input_shape": {"class_name": "TensorShape", "items": [null, 1]}}}, "training_config": {"loss": "mse", "metrics": ["mse"], "weighted_metrics": null, "loss_weights": null, "sample_weight_mode": null, "optimizer_config": {"class_name": "Adadelta", "config": {"name": "Adadelta", "learning_rate": 0.0010000000474974513, "decay": 0.0, "rho": 0.949999988079071, "epsilon": 1e-07}}}}
�

kernel
bias
regularization_losses
	variables
trainable_variables
	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"class_name": "Dense", "name": "dense_874", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": {"class_name": "__tuple__", "items": [null, 1]}, "stateful": false, "config": {"name": "dense_874", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 1]}, "dtype": "float32", "units": 40, "activation": "linear", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 1}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 1]}}
�

kernel
bias
regularization_losses
	variables
trainable_variables
	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"class_name": "Dense", "name": "dense_875", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "config": {"name": "dense_875", "trainable": true, "dtype": "float32", "units": 40, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 40}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 40]}}
�

kernel
bias
regularization_losses
	variables
trainable_variables
	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"class_name": "Dense", "name": "dense_876", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "config": {"name": "dense_876", "trainable": true, "dtype": "float32", "units": 40, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 40}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 40]}}
�

 kernel
!bias
"regularization_losses
#	variables
$trainable_variables
%	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"class_name": "Dense", "name": "dense_877", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "config": {"name": "dense_877", "trainable": true, "dtype": "float32", "units": 40, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 40}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 40]}}
�

&kernel
'bias
(regularization_losses
)	variables
*trainable_variables
+	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"class_name": "Dense", "name": "dense_878", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "config": {"name": "dense_878", "trainable": true, "dtype": "float32", "units": 40, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 40}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 40]}}
�

,kernel
-bias
.regularization_losses
/	variables
0trainable_variables
1	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"class_name": "Dense", "name": "dense_879", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "config": {"name": "dense_879", "trainable": true, "dtype": "float32", "units": 40, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 40}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 40]}}
�

2kernel
3bias
4regularization_losses
5	variables
6trainable_variables
7	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"class_name": "Dense", "name": "dense_880", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "config": {"name": "dense_880", "trainable": true, "dtype": "float32", "units": 1, "activation": "linear", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 40}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 40]}}
�
8iter
	9decay
:learning_rate
;rho
accum_grado
accum_gradp
accum_gradq
accum_gradr
accum_grads
accum_gradt 
accum_gradu!
accum_gradv&
accum_gradw'
accum_gradx,
accum_grady-
accum_gradz2
accum_grad{3
accum_grad|	accum_var}	accum_var~	accum_var	accum_var�	accum_var�	accum_var� 	accum_var�!	accum_var�&	accum_var�'	accum_var�,	accum_var�-	accum_var�2	accum_var�3	accum_var�"
	optimizer
 "
trackable_list_wrapper
�
0
1
2
3
4
5
 6
!7
&8
'9
,10
-11
212
313"
trackable_list_wrapper
�
0
1
2
3
4
5
 6
!7
&8
'9
,10
-11
212
313"
trackable_list_wrapper
�
<layer_metrics
	regularization_losses
=layer_regularization_losses
>non_trainable_variables

?layers
@metrics

	variables
trainable_variables
�__call__
�_default_save_signature
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
-
�serving_default"
signature_map
": (2dense_874/kernel
:(2dense_874/bias
 "
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
�
Alayer_metrics
regularization_losses
Blayer_regularization_losses
Cnon_trainable_variables

Dlayers
Emetrics
	variables
trainable_variables
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
": ((2dense_875/kernel
:(2dense_875/bias
 "
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
�
Flayer_metrics
regularization_losses
Glayer_regularization_losses
Hnon_trainable_variables

Ilayers
Jmetrics
	variables
trainable_variables
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
": ((2dense_876/kernel
:(2dense_876/bias
 "
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
�
Klayer_metrics
regularization_losses
Llayer_regularization_losses
Mnon_trainable_variables

Nlayers
Ometrics
	variables
trainable_variables
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
": ((2dense_877/kernel
:(2dense_877/bias
 "
trackable_list_wrapper
.
 0
!1"
trackable_list_wrapper
.
 0
!1"
trackable_list_wrapper
�
Player_metrics
"regularization_losses
Qlayer_regularization_losses
Rnon_trainable_variables

Slayers
Tmetrics
#	variables
$trainable_variables
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
": ((2dense_878/kernel
:(2dense_878/bias
 "
trackable_list_wrapper
.
&0
'1"
trackable_list_wrapper
.
&0
'1"
trackable_list_wrapper
�
Ulayer_metrics
(regularization_losses
Vlayer_regularization_losses
Wnon_trainable_variables

Xlayers
Ymetrics
)	variables
*trainable_variables
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
": ((2dense_879/kernel
:(2dense_879/bias
 "
trackable_list_wrapper
.
,0
-1"
trackable_list_wrapper
.
,0
-1"
trackable_list_wrapper
�
Zlayer_metrics
.regularization_losses
[layer_regularization_losses
\non_trainable_variables

]layers
^metrics
/	variables
0trainable_variables
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
": (2dense_880/kernel
:2dense_880/bias
 "
trackable_list_wrapper
.
20
31"
trackable_list_wrapper
.
20
31"
trackable_list_wrapper
�
_layer_metrics
4regularization_losses
`layer_regularization_losses
anon_trainable_variables

blayers
cmetrics
5	variables
6trainable_variables
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
:	 (2Adadelta/iter
: (2Adadelta/decay
 : (2Adadelta/learning_rate
: (2Adadelta/rho
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
Q
0
1
2
3
4
5
6"
trackable_list_wrapper
.
d0
e1"
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
�
	ftotal
	gcount
h	variables
i	keras_api"�
_tf_keras_metricj{"class_name": "Mean", "name": "loss", "dtype": "float32", "config": {"name": "loss", "dtype": "float32"}}
�
	jtotal
	kcount
l
_fn_kwargs
m	variables
n	keras_api"�
_tf_keras_metric�{"class_name": "MeanMetricWrapper", "name": "mse", "dtype": "float32", "config": {"name": "mse", "dtype": "float32", "fn": "mean_squared_error"}}
:  (2total
:  (2count
.
f0
g1"
trackable_list_wrapper
-
h	variables"
_generic_user_object
:  (2total
:  (2count
 "
trackable_dict_wrapper
.
j0
k1"
trackable_list_wrapper
-
m	variables"
_generic_user_object
4:2(2$Adadelta/dense_874/kernel/accum_grad
.:,(2"Adadelta/dense_874/bias/accum_grad
4:2((2$Adadelta/dense_875/kernel/accum_grad
.:,(2"Adadelta/dense_875/bias/accum_grad
4:2((2$Adadelta/dense_876/kernel/accum_grad
.:,(2"Adadelta/dense_876/bias/accum_grad
4:2((2$Adadelta/dense_877/kernel/accum_grad
.:,(2"Adadelta/dense_877/bias/accum_grad
4:2((2$Adadelta/dense_878/kernel/accum_grad
.:,(2"Adadelta/dense_878/bias/accum_grad
4:2((2$Adadelta/dense_879/kernel/accum_grad
.:,(2"Adadelta/dense_879/bias/accum_grad
4:2(2$Adadelta/dense_880/kernel/accum_grad
.:,2"Adadelta/dense_880/bias/accum_grad
3:1(2#Adadelta/dense_874/kernel/accum_var
-:+(2!Adadelta/dense_874/bias/accum_var
3:1((2#Adadelta/dense_875/kernel/accum_var
-:+(2!Adadelta/dense_875/bias/accum_var
3:1((2#Adadelta/dense_876/kernel/accum_var
-:+(2!Adadelta/dense_876/bias/accum_var
3:1((2#Adadelta/dense_877/kernel/accum_var
-:+(2!Adadelta/dense_877/bias/accum_var
3:1((2#Adadelta/dense_878/kernel/accum_var
-:+(2!Adadelta/dense_878/bias/accum_var
3:1((2#Adadelta/dense_879/kernel/accum_var
-:+(2!Adadelta/dense_879/bias/accum_var
3:1(2#Adadelta/dense_880/kernel/accum_var
-:+2!Adadelta/dense_880/bias/accum_var
�2�
/__inference_sequential_87_layer_call_fn_1017132
/__inference_sequential_87_layer_call_fn_1016851
/__inference_sequential_87_layer_call_fn_1017099
/__inference_sequential_87_layer_call_fn_1016923�
���
FullArgSpec1
args)�&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�2�
J__inference_sequential_87_layer_call_and_return_conditional_losses_1017066
J__inference_sequential_87_layer_call_and_return_conditional_losses_1017015
J__inference_sequential_87_layer_call_and_return_conditional_losses_1016778
J__inference_sequential_87_layer_call_and_return_conditional_losses_1016739�
���
FullArgSpec1
args)�&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�2�
"__inference__wrapped_model_1016547�
���
FullArgSpec
args� 
varargsjargs
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *.�+
)�&
dense_874_input���������
�2�
+__inference_dense_874_layer_call_fn_1017151�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
F__inference_dense_874_layer_call_and_return_conditional_losses_1017142�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
+__inference_dense_875_layer_call_fn_1017171�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
F__inference_dense_875_layer_call_and_return_conditional_losses_1017162�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
+__inference_dense_876_layer_call_fn_1017191�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
F__inference_dense_876_layer_call_and_return_conditional_losses_1017182�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
+__inference_dense_877_layer_call_fn_1017211�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
F__inference_dense_877_layer_call_and_return_conditional_losses_1017202�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
+__inference_dense_878_layer_call_fn_1017231�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
F__inference_dense_878_layer_call_and_return_conditional_losses_1017222�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
+__inference_dense_879_layer_call_fn_1017251�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
F__inference_dense_879_layer_call_and_return_conditional_losses_1017242�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
+__inference_dense_880_layer_call_fn_1017270�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
F__inference_dense_880_layer_call_and_return_conditional_losses_1017261�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
<B:
%__inference_signature_wrapper_1016964dense_874_input�
"__inference__wrapped_model_1016547� !&',-238�5
.�+
)�&
dense_874_input���������
� "5�2
0
	dense_880#� 
	dense_880����������
F__inference_dense_874_layer_call_and_return_conditional_losses_1017142\/�,
%�"
 �
inputs���������
� "%�"
�
0���������(
� ~
+__inference_dense_874_layer_call_fn_1017151O/�,
%�"
 �
inputs���������
� "����������(�
F__inference_dense_875_layer_call_and_return_conditional_losses_1017162\/�,
%�"
 �
inputs���������(
� "%�"
�
0���������(
� ~
+__inference_dense_875_layer_call_fn_1017171O/�,
%�"
 �
inputs���������(
� "����������(�
F__inference_dense_876_layer_call_and_return_conditional_losses_1017182\/�,
%�"
 �
inputs���������(
� "%�"
�
0���������(
� ~
+__inference_dense_876_layer_call_fn_1017191O/�,
%�"
 �
inputs���������(
� "����������(�
F__inference_dense_877_layer_call_and_return_conditional_losses_1017202\ !/�,
%�"
 �
inputs���������(
� "%�"
�
0���������(
� ~
+__inference_dense_877_layer_call_fn_1017211O !/�,
%�"
 �
inputs���������(
� "����������(�
F__inference_dense_878_layer_call_and_return_conditional_losses_1017222\&'/�,
%�"
 �
inputs���������(
� "%�"
�
0���������(
� ~
+__inference_dense_878_layer_call_fn_1017231O&'/�,
%�"
 �
inputs���������(
� "����������(�
F__inference_dense_879_layer_call_and_return_conditional_losses_1017242\,-/�,
%�"
 �
inputs���������(
� "%�"
�
0���������(
� ~
+__inference_dense_879_layer_call_fn_1017251O,-/�,
%�"
 �
inputs���������(
� "����������(�
F__inference_dense_880_layer_call_and_return_conditional_losses_1017261\23/�,
%�"
 �
inputs���������(
� "%�"
�
0���������
� ~
+__inference_dense_880_layer_call_fn_1017270O23/�,
%�"
 �
inputs���������(
� "�����������
J__inference_sequential_87_layer_call_and_return_conditional_losses_1016739y !&',-23@�=
6�3
)�&
dense_874_input���������
p

 
� "%�"
�
0���������
� �
J__inference_sequential_87_layer_call_and_return_conditional_losses_1016778y !&',-23@�=
6�3
)�&
dense_874_input���������
p 

 
� "%�"
�
0���������
� �
J__inference_sequential_87_layer_call_and_return_conditional_losses_1017015p !&',-237�4
-�*
 �
inputs���������
p

 
� "%�"
�
0���������
� �
J__inference_sequential_87_layer_call_and_return_conditional_losses_1017066p !&',-237�4
-�*
 �
inputs���������
p 

 
� "%�"
�
0���������
� �
/__inference_sequential_87_layer_call_fn_1016851l !&',-23@�=
6�3
)�&
dense_874_input���������
p

 
� "�����������
/__inference_sequential_87_layer_call_fn_1016923l !&',-23@�=
6�3
)�&
dense_874_input���������
p 

 
� "�����������
/__inference_sequential_87_layer_call_fn_1017099c !&',-237�4
-�*
 �
inputs���������
p

 
� "�����������
/__inference_sequential_87_layer_call_fn_1017132c !&',-237�4
-�*
 �
inputs���������
p 

 
� "�����������
%__inference_signature_wrapper_1016964� !&',-23K�H
� 
A�>
<
dense_874_input)�&
dense_874_input���������"5�2
0
	dense_880#� 
	dense_880���������