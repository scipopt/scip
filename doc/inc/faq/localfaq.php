<!-- this file is used to generate the local doxygen documentation -->
<!-- using make doc from within scip or soplex -->

<style>
.reveal:hover {
    text-shadow: 1px 1px 1px #777;
}
.answer {
    background-color: #fff;
    padding-left:   1em;
}
</style>

<?php include('faqdata.php'); ?>
<?php
//Output a table of contents with popping up

$sectionCounter = 1;
foreach ($faq as $section) {
  echo '<h3>'.$section['title'].'</h3>';
  echo '<ol>';
    foreach($section['content'] as $item) {
    $label = $item['label'];
  ?>
  <li>
    <div id="<?php echo $label ?>" class="targetpadding">
      <div class="reveal_faq">
        <a href="#<?php echo $label ?>">
            <?php echo '<h4>'.$item['question'].'</h4>' ?>
        </a>
        </div>
    </div>
    <div id="<?php echo $label ?>_ans" class="answer">
       <?php echo $item['answer'];?>
    </div>
  </li>
  <?php
  }
  echo "</ol><br/>\n";
}
?>
<!--WE don't use this stuff because it is uncontrollable
<script>
$(".reveal").click(function() {
    var questionId = "#"+this.id+"_ans";
    $(questionId).toggle("fast");
});
</script>-->
