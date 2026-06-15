# Override miniconda Python with container's Python
export PATH=/opt/spack/opt/spack/__spack_path_placeholder__/__spack_path_placeholder__/__spack_path_placeholder__/__spack_path_placeholder__/linux-x86_64/python-3.12.9-rjkgjb7l4xizmkw2ivgug3gajs2qimux/bin:$PATH

# Custom LCContent
export MARLIN_DLL=$HOME/myLCContent/lib64/libLCContent.so:$MARLIN_DLL

# Library paths
export LD_LIBRARY_PATH=/opt/spack/opt/spack/__spack_path_placeholder__/__spack_path_placeholder__/__spack_path_placeholder__/__spack_path_placeholder__/linux-x86_64/pandorasdk-3.4.2-gdgqi3cx2lkvqn3mqtztyarhkoxg7hm6/lib:/opt/spack/opt/spack/__spack_path_placeholder__/__spack_path_placeholder__/__spack_path_placeholder__/__spack_path_placeholder__/linux-x86_64/pandoramonitoring-3.6.0-bxk4em6f2apbjcva7lghcllzcuny4oqj/lib:$(root-config --libdir):$LD_LIBRARY_PATH
