#include "RecoE2E/FrameProducers/interface/predict_tf.h"
#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"
std::string print_shape(const std::vector<int64_t>& v) {
  std::stringstream ss("");
  for (size_t i = 0; i < v.size() - 1; i++)
    ss << v[i] << "x";
  ss << v[v.size() - 1];
  return ss.str();
}
template <typename T>
T vectorProduct(const std::vector<T>& v)
{
    return accumulate(v.begin(), v.end(), 1, std::multiplies<T>());
}

e2e::Frame2D e2e::predict_tf(e2e::Frame4D& vinputFrame, string model_filename, string input_layer_name, string output_layer_name){
 //std::cout<<"In onnx predict.cc";
 e2e::Frame2D output_preds;
 tensorflow::Session* session;
 tensorflow::GraphDef graph_def;
 tensorflow::SessionOptions opts;
 std::vector<tensorflow::Tensor> outputs; // Store outputs

 //Onxx
 std::unique_ptr<cms::Ort::ONNXRuntime> model;
 Ort::SessionOptions session_options;
 std::vector<Ort::Value> inputTensors;
 std::vector<Ort::Value> outputTensors;
 auto providers = Ort::GetAvailableProviders();
 
 //session_options.AppendExecutionProvider_CUDA(cuda_options);
 // Sets graph optimization level
 session_options.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_EXTENDED);
 std::string instanceName{"Particle-classification-inference"};
 Ort::Env env(OrtLoggingLevel::ORT_LOGGING_LEVEL_WARNING,
                 instanceName.c_str());
 std::string modelFilepath ="RecoE2E/"+model_filename;
 OrtCUDAProviderOptions cuda_options{};
 session_options.AppendExecutionProvider_CUDA(cuda_options);
 Ort::Session ort_session(env, modelFilepath.c_str(), session_options);

// for (auto provider : providers) {
//     std::cout << "Onxx Provider : "<<provider << endl;
//    // std::cout<<session_options.config.gpu_options().visible_device_list();
// }
const int64_t batchSize = 2;

Ort::AllocatorWithDefaultOptions allocator;

    size_t numInputNodes = ort_session.GetInputCount();
    size_t numOutputNodes = ort_session.GetOutputCount();

    const char* inputName = ort_session.GetInputName(0, allocator);

    Ort::TypeInfo inputTypeInfo = ort_session.GetInputTypeInfo(0);
    auto inputTensorInfo = inputTypeInfo.GetTensorTypeAndShapeInfo();

    ONNXTensorElementDataType inputType = inputTensorInfo.GetElementType();

    std::vector<int64_t> inputDims = inputTensorInfo.GetShape();
    if (inputDims.at(0) == -1)
    {
        // std::cout << "Got dynamic batch size. Setting input batch size to "
        //           << vinputFrame.size() << "." << std::endl;
        inputDims.at(0) = vinputFrame.size();
    }

    const char* outputName = ort_session.GetOutputName(0, allocator);

    Ort::TypeInfo outputTypeInfo = ort_session.GetOutputTypeInfo(0);
    auto outputTensorInfo = outputTypeInfo.GetTensorTypeAndShapeInfo();

    ONNXTensorElementDataType outputType = outputTensorInfo.GetElementType();

    std::vector<int64_t> outputDims = outputTensorInfo.GetShape();
    if (outputDims.at(0) == -1)
    {
        // std::cout << "Got dynamic batch size. Setting output batch size to "
        //           << vinputFrame.size() << "." << std::endl;
        outputDims.at(0) = vinputFrame.size();
    }

    // std::cout << "Number of Input Nodes: " << numInputNodes << std::endl;
    // std::cout << "Number of Output Nodes: " << numOutputNodes << std::endl;
    // std::cout << "Input Name: " << inputName << std::endl;
    // std::cout << "Input Type: " << inputType << std::endl;
    // std::cout << "Input Dimensions: " << print_shape(inputDims) << std::endl;
    // std::cout << "Output Name: " << outputName << std::endl;
    // std::cout << "Output Type: " << outputType << std::endl;
    //  std::cout << "Output Dimensions: " << print_shape(outputDims) << std::endl;
    auto inputShape = ort_session.GetInputTypeInfo(0).GetTensorTypeAndShapeInfo().GetShape();

    int batch_sz = vinputFrame.size();
    int no_channels = vinputFrame[0].size();
    int frame_width = vinputFrame[0][0][0].size();
    int frame_height = vinputFrame[0][0].size();
    //const std::array<int64_t, 4> inputTensorValues = { batch_sz, no_channels, frame_width, frame_height };
    // std::vector inputTensorValues(batch_sz, vector<vector<vector<float>>>
    //                                        (no_channels, vector<vector<float>>
    //                                        (frame_width, vector<float>
    //                                        (frame_height))));
    //std::vector inputTensorValues(batch_sz*no_channels*frame_width*frame_height, 0);
    size_t inputTensorSize = vectorProduct(inputDims);
    std::vector<float> inputTensorValues(inputTensorSize);
    //std::cout << "Vector product : " << inputTensorSize << std::endl;
    //float inputTensorValues[batch_sz][no_channels][frame_width][frame_height];
    int c =0;
    for (int batch_idx=0; batch_idx<batch_sz; batch_idx++){
      for (int row_idx=0;row_idx<frame_height; row_idx++){
       for (int col_idx=0; col_idx<frame_width; col_idx++){
        for (int depth_idx=0; depth_idx<no_channels; depth_idx++){
         inputTensorValues[c] = vinputFrame[batch_idx][depth_idx][row_idx][col_idx];
         c++;
        }
       }
      }
     }
     std::vector<const char*> inputNames{inputName};
    std::vector<const char*> outputNames{outputName};
     size_t outputTensorSize = vectorProduct(outputDims);
     std::vector<float> outputTensorValues(outputTensorSize);
     
    // define Tensor
    auto memory_info = Ort::MemoryInfo::CreateCpu(OrtDeviceAllocator, OrtMemTypeCPU);
    //auto inputTensor = Ort::Value::CreateTensor<float>(memory_info, inputTensorValues.data(), inputTensorValues.size(), inputShape.data(), inputShape.size());
    //std::vector<Ort::Value> inputTensors;
    inputTensors.push_back(Ort::Value::CreateTensor<float>(
        memory_info, inputTensorValues.data(), inputTensorSize, inputDims.data(),
        inputDims.size()));

    outputTensors.push_back(Ort::Value::CreateTensor<float>(
        memory_info, outputTensorValues.data(), outputTensorSize,
        outputDims.data(), outputDims.size()));

    ort_session.Run(Ort::RunOptions{nullptr}, inputNames.data(),
                inputTensors.data(), 1 /*Number of inputs*/, outputNames.data(),
                outputTensors.data(), 1 /*Number of outputs*/);


    //input_tensors.push_back(Ort::Value::CreateTensor<float>(memory_info,inputTensorValues.data(), inputTensorValues.size(), inputShape.data(), inputShape.size()));
    //auto outputTensor = Ort::Value::CreateTensor<float>(memory_info, results.data(), results.size(), outputShape.data(), outputShape.size());
    //const std::array<int64_t, 2> outputShape = { 1, numClasses };
    //std::vector<vector<vector<vector<float>>>> inputTensorValues(batch_sz, vector<float> (no_channels, vector<float> (frame_width, vector<float> (frame_height, 0))));

    //auto allocator_info = Ort::AllocatorInfo::CreateCpu(OrtDeviceAllocator, OrtMemTypeCPU);

    //std::cout<<"Onxx Over"<<std::endl;
    //std::cout<<outputTensorValues.at(0)<<std::endl;
//  create a new session
//  TF_CHECK_OK(NewSession(opts, &session));
 
//  std::string graph_definition="RecoE2E/"+model_filename;
 //std::cout<<" >> Running Inference. 1 "<<std::endl;
 edm::LogInfo("JetFrameProducer") << " >> Running Inference.";
 //int batch_sz = vinputFrame.size();
 if (batch_sz>0){ 
  //int no_channels = vinputFrame[0].size();
  //std::cout<<no_channels<<std::endl;
  edm::LogInfo("JetFrameProducer") << no_channels;
  if (no_channels>0){
   //int frame_height = vinputFrame[0][0].size();
   if (frame_height>0){
    //int frame_width = vinputFrame[0][0][0].size();
    if (frame_width>0){
       // std::cout<<" >> Running Inference."<<std::endl;
     //TF_CHECK_OK(ReadBinaryProto(Env::Default(), graph_definition, &graph_def));
     // load the graph definition, i.e. an object that contains the computational graph
    //  tensorflow::GraphDef* graphDef = tensorflow::loadGraphDef(graph_definition);
    //  tensorflow::Tensor tmp(tensorflow::DT_FLOAT, tensorflow::TensorShape({frame_height, frame_width}));
  
    //  tensorflow::Tensor x(tensorflow::DT_FLOAT, tensorflow::TensorShape({batch_sz, frame_height, frame_width,  no_channels}));
    //  auto _XTensor = x.tensor<float,4>();
    //  for (int batch_idx=0; batch_idx<batch_sz; batch_idx++){
    //   for (int row_idx=0;row_idx<frame_height; row_idx++){
    //    for (int col_idx=0; col_idx<frame_width; col_idx++){
    //     for (int depth_idx=0; depth_idx<no_channels; depth_idx++){
    //      _XTensor(batch_idx, row_idx, col_idx, depth_idx) = vinputFrame[batch_idx][depth_idx][row_idx][col_idx];
    //     }
    //    }
    //   }
    //  }
     /*if(!x.CopyFrom(tmp, tensorflow::TensorShape({1, frame_height, frame_width, 1}))){
      std::cout<<" >> Reshape not successfull."<<endl;
     }*/
     // Set GPU options
    //  tensorflow::graph::SetDefaultDevice("/gpu:0", &graph_def);
    //  opts.config.mutable_gpu_options()->set_per_process_gpu_memory_fraction(0.5);
    //  opts.config.mutable_gpu_options()->set_allow_growth(true);
 
    //  //int GPUID = std::stoi(params->getGpuDeviceStr());
    //  int GPUID = 0;
    //  setenv("CUDA_VISIBLE_DEVICES", "", GPUID);

    //  //std::cout << "Initial  visible_device_list : "<<opts.config.gpu_options().visible_device_list() << std::endl;
    //  opts.config.mutable_gpu_options()->set_allow_growth(true);
    //  opts.config.mutable_gpu_options()->set_per_process_gpu_memory_fraction(0.5);//params->getGpuMemoryRatio());
 
 
    //  // Load graph into session
    //  //TF_CHECK_OK(session->Create(graph_def));
 
    //  // create a session
    //  session = tensorflow::createSession(graphDef);
 
     // Initialize our variables
 
     //TF_CHECK_OK(session->Run({{input_layer_name/*"inputs"*/, x}/*, {"y", y}*/}, {output_layer_name/*"softmax_1/Sigmoid"*/}, {}, &outputs)); // Get output
     //tensorflow::run(session, { { "x", x }, {"y", y} }, { "cost" }, &outputs);
     //std::cout<<" >> Classification predictions: "<<std::endl;
     edm::LogInfo("JetFrameProducer") << " >> Classification predictions: ";     
    //std::cout<<"Batch size : "<<batch_sz<<std::endl;
     // if (outputTensorValues[0].shape().dims()!=2) edm::LogInfo("JetFrameProducer") << " * Expected 2 dimensional output. Received " << outputs[0].shape().dims() << " dimensional output."; //std::cout<<"* Expected 2 dimensional output. Received "<<outputs[0].shape().dims()<<" dimension output."<<std::endl;
    if(batch_sz==0) edm::LogInfo("JetFrameProducer") << " * Expected 2";
     else {
        //std::cout<<"before serializing outputs "<<std::endl;
        int co =0;
      e2e::Frame1D preds;
      for (int row_idx=0; row_idx<batch_sz; row_idx++){
       for (int col_idx=0; col_idx<1; col_idx++){
        //std::cout<<"Inside for loop 2 "<<std::endl;
        preds.push_back(outputTensorValues.at(co));
        co++;
       }
       //std::cout<<"Inside for loop 1 "<<std::endl;
       output_preds.push_back(preds);
      }
      //std::cout<<"Push back output "<<std::endl;
      for (int row_idx=0; row_idx<batch_sz; row_idx++){
       for (int col_idx=0; col_idx<1; col_idx++){
        //std::cout<<"outputs: ("<<row_idx<<","<<col_idx<<"): "<<output_preds[row_idx][col_idx]<<std::endl;
	edm::LogInfo("JetFrameProducer") << "outputs: (" << row_idx << "," << col_idx << "): " << output_preds[row_idx][col_idx];
       }
      }
     }
     //std::cout<<" >> Batch size is: "<<outputs[0].shape().dim_size(0)<<std::endl;
    //  edm::LogInfo("JetFrameProducer") << " >> Batch size is: " << outputs[0].shape().dim_size(0);
    //  outputs.clear();
  
    //  session->Close();
    //  delete session;
     //std::cout<<" >> Classification done"<<endl;
     edm::LogInfo("JetFrameProducer") << " >> Classification done.";
    }
    else{
     //std::cout<<"* Shape Error: Invalid Width(<0) dimension. Expected format: (N, C, H, W)"<<std::endl;
     edm::LogInfo("JetFrameProducer") << " * Shape Error: Invalid Width (<0) dimension. Expected format: (N, C, H, W)";	
    }
   }
   else{
    //std::cout<<"* Shape Error: Invalid Height(<0) dimension. Expected format: (N, C, H, W)"<<std::endl;
    edm::LogInfo("JetFrameProducer") << " * Shape Error: Invalid Height (<0) dimension. Expected format: (N, C, H, W)";
   }
  }
  else{
   //std::cout<<"* Shape Error: Invalid Channel(<0) dimension. Expected format: (N, C, H, W)"<<std::endl;
   edm::LogInfo("JetFrameProducer") << " * Shape Error: Invalid Channel (<0) dimension. Expected format: (N, C, H, W)";
  }
 }
 else{
  //std::cout<<"* Shape Error: Invalid Batch(<0) dimension. Expected format: (N, C, H, W)"<<std::endl;
  edm::LogInfo("JetFrameProducer") << " * Shape Error: Invalid Batch (<0) dimension. Expected format: (N, C, H, W)";
 }
 // cleanup
 //tensorflow::closeSession(session);
 //delete graphDef;
 return output_preds;
}
