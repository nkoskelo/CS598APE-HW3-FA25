docker run -it   --cap-add SYS_ADMIN \
  --cap-add SYS_PTRACE \
  --cap-add PERFMON --security-opt seccomp=unconfined -v  `pwd`:/host wsmoses/598ape /bin/bash

